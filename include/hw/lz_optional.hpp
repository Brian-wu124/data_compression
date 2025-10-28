/*
 * (c) Copyright 2019-2022 Xilinx, Inc. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#ifndef _XFCOMPRESSION_LZ_OPTIONAL_HPP_
#define _XFCOMPRESSION_LZ_OPTIONAL_HPP_

/**
 * @file lz_optional.hpp
 * @brief Header for modules used in LZ4 and snappy compression kernels.
 *
 * This file is part of Vitis Data Compression Library.
 */

#include "compress_utils.hpp"
#include "hls_stream.h"

#include <ap_int.h>
#include <assert.h>
#include <stdint.h>

// Dynamic Huffman Frequency counts
// Based on GZIP/ZLIB spec
#include "zlib_tables.hpp"

#define d_code(dist, dist_code) ((dist) < 256 ? dist_code[dist] : dist_code[256 + ((dist) >> 7)])

namespace xf {
namespace compression {

typedef ap_uint<32> compressd_dt;

static void lz77Divide(hls::stream<compressd_dt>& inStream,
                       hls::stream<ap_uint<32> >& outStream,
                       hls::stream<bool>& endOfStream,
                       hls::stream<uint32_t>& outStreamTree,
                       hls::stream<uint32_t>& compressedSize,
                       uint32_t input_size) {
    if (input_size == 0) {
        compressedSize << 0;
        outStream << 0;
        endOfStream << 1;
        outStreamTree << 9999;
        return;
    }

    const uint16_t c_ltree_size = 286;
    const uint16_t c_dtree_size = 30;

    uint32_t lcl_dyn_ltree[c_ltree_size];
    uint32_t lcl_dyn_dtree[c_dtree_size];

ltree_init:
    for (uint32_t i = 0; i < c_ltree_size; i++) {
#pragma HLS PIPELINE off
        lcl_dyn_ltree[i] = 0;
    }
dtree_init:
    for (uint32_t i = 0; i < c_dtree_size; i++) {
#pragma HLS PIPELINE off
        lcl_dyn_dtree[i] = 0;
    }
    int length = 0;
    int code = 0;
    int n = 0;

    uint32_t out_cntr = 0;
    uint8_t match_len = 0;
    uint32_t loc_idx = 0;

    compressd_dt nextEncodedValue = inStream.read();
    uint32_t cSizeCntr = 0;
lz77_divide:
    for (uint32_t i = 0; i < input_size; i++) {
#pragma HLS LOOP_TRIPCOUNT min = 1048576 max = 1048576
#pragma HLS PIPELINE II = 1
#pragma HLS dependence variable = lcl_dyn_ltree inter false
#pragma HLS dependence variable = lcl_dyn_dtree inter false
        compressd_dt tmpEncodedValue = nextEncodedValue;
        if (i < (input_size - 1)) nextEncodedValue = inStream.read();

        uint8_t tCh = tmpEncodedValue.range(7, 0);
        uint8_t tLen = tmpEncodedValue.range(15, 8);
        uint16_t tOffset = tmpEncodedValue.range(31, 16);
        uint32_t ltreeIdx, dtreeIdx;
        if (tLen > 0) {
            ltreeIdx = length_code[tLen] + 256 + 1;
            dtreeIdx = d_code(tOffset, dist_code);

            i += (tLen - 1);
            tmpEncodedValue.range(15, 8) = tLen;
            lcl_dyn_ltree[ltreeIdx]++;
            lcl_dyn_dtree[dtreeIdx]++;
        } else {
            ltreeIdx = tCh;
            lcl_dyn_ltree[ltreeIdx]++;
        }

        outStream << tmpEncodedValue;
        endOfStream << 0;
        cSizeCntr++;
    }

    compressedSize << (cSizeCntr * 4);
    outStream << 0;
    endOfStream << 1;

    for (uint32_t i = 0; i < c_ltree_size; i++) {
#pragma HLS PIPELINE off
        outStreamTree << lcl_dyn_ltree[i];
    }
    for (uint32_t i = 0; i < c_dtree_size; i++) {
#pragma HLS PIPELINE off
        outStreamTree << lcl_dyn_dtree[i];
    }
}

template <int MAX_FREQ_DWIDTH = 32, int CORE_ID = 0>
static void lz77DivideStream(hls::stream<IntVectorStream_dt<32, 1> >& inStream,
                             hls::stream<IntVectorStream_dt<9, 1> >& outStream,
                             hls::stream<IntVectorStream_dt<MAX_FREQ_DWIDTH, 1> >& outTreeStream) {
    constexpr uint16_t c_ltree_size = 286;
    constexpr uint16_t c_dtree_size = 30;
    // lz77 encoder states
    enum LZ77EncoderStates { WRITE_LITERAL, WRITE_OFFSET0, WRITE_OFFSET1 };
    IntVectorStream_dt<9, 1> outValue;
    IntVectorStream_dt<MAX_FREQ_DWIDTH, 1> outTreeVal;
    bool last_block = false;
    bool just_started = true;
#ifndef AVOID_STATIC_MODE
    static bool resetTree = true;
#else
    bool resetTree = true;
#endif
    ap_uint<MAX_FREQ_DWIDTH> lcl_dyn_ltree[c_ltree_size];
    ap_uint<MAX_FREQ_DWIDTH> lcl_dyn_dtree[c_dtree_size];
#pragma HLS bind_storage variable = lcl_dyn_ltree type = ram_2p impl = LUTRAM
#pragma HLS bind_storage variable = lcl_dyn_dtree type = ram_2p impl = LUTRAM
ltree_init:
    for (uint16_t i = 0; i < c_ltree_size && resetTree; i++) lcl_dyn_ltree[i] = 0;
dtree_init:
    for (uint16_t i = 0; i < c_dtree_size && resetTree; i++) lcl_dyn_dtree[i] = 0;
    resetTree = false;

    while (!last_block) { // iterate over multiple blocks in a file
        enum LZ77EncoderStates next_state = WRITE_LITERAL;

        uint8_t tCh = 0;
        uint8_t tLen = 0;
        uint8_t tOffset1 = 0;
        uint8_t tOffset2 = 0;
        // set ouput data to be 1 valid length
        outValue.strobe = 1;
    lz77_divide:
        while (true) {
#pragma HLS LOOP_TRIPCOUNT min = 1048576 max = 1048576
#pragma HLS PIPELINE II = 1
#pragma HLS dependence variable = lcl_dyn_ltree inter false
#pragma HLS dependence variable = lcl_dyn_dtree inter false
            uint32_t ltreeIdx, dtreeIdx;
            if (next_state == WRITE_LITERAL) {
                auto encodedValue = inStream.read();
                if (encodedValue.strobe == 0) {
                    last_block = just_started;
                    just_started = true;
                    break;
                }
                just_started = false;
                tCh = encodedValue.data[0].range(7, 0);
                tLen = encodedValue.data[0].range(15, 8);
                uint16_t tOffset = encodedValue.data[0].range(31, 16);
                tOffset1 = encodedValue.data[0].range(23, 16);
                tOffset2 = encodedValue.data[0].range(31, 24);

                if (tLen) {
                    ltreeIdx = length_code[tLen] + 256 + 1;
                    dtreeIdx = d_code(tOffset, dist_code);
                    lcl_dyn_ltree[ltreeIdx]++;
                    lcl_dyn_dtree[dtreeIdx]++;
                    outValue.data[0].range(8, 8) = 1;
                    outValue.data[0].range(7, 0) = tLen;
                    next_state = WRITE_OFFSET0;
                } else {
                    ltreeIdx = tCh;
                    lcl_dyn_ltree[ltreeIdx]++;
                    outValue.data[0] = tCh;
                    next_state = WRITE_LITERAL;
                }
            } else if (next_state == WRITE_OFFSET0) {
                outValue.data[0] = tOffset1;
                next_state = WRITE_OFFSET1;
            } else if (next_state == WRITE_OFFSET1) {
                outValue.data[0] = tOffset2;
                next_state = WRITE_LITERAL;
            }
            outStream << outValue;
        }
        // File/Block end Condition
        outValue.data[0] = 0;
        outValue.strobe = 0;
        outStream << outValue;
        // write literal and distance trees
        if (!last_block) {
            outTreeVal.strobe = 1;
        read_write_ltree:
            for (ap_uint<9> i = 0; i < c_ltree_size; ++i) {
#pragma HLS PIPELINE II = 1
                outTreeVal.data[0] = lcl_dyn_ltree[i];
                lcl_dyn_ltree[i] = 0;
                outTreeStream << outTreeVal;
            }
        read_write_dtree:
            for (ap_uint<9> i = 0; i < c_dtree_size; ++i) {
#pragma HLS PIPELINE II = 1
                outTreeVal.data[0] = lcl_dyn_dtree[i];
                lcl_dyn_dtree[i] = 0;
                outTreeStream << outTreeVal;
            }
        }
    }
    // eos needed only to indicated end of block
    outTreeVal.strobe = 0;
    outTreeStream << outTreeVal;
}

static void lz77DivideStatic(hls::stream<IntVectorStream_dt<32, 1> >& inStream,
                             hls::stream<IntVectorStream_dt<9, 1> >& outStream) {
    // lz77 encoder states
    enum LZ77EncoderStates { WRITE_LITERAL, WRITE_OFFSET0, WRITE_OFFSET1 };
    IntVectorStream_dt<9, 1> outValue;
    bool last_block = false;
    bool just_started = true;

    while (!last_block) { // iterate over multiple blocks in a file
        enum LZ77EncoderStates next_state = WRITE_LITERAL;

        uint8_t tCh = 0;
        uint8_t tLen = 0;
        uint8_t tOffset1 = 0;
        uint8_t tOffset2 = 0;
        // set output data to be valid
        outValue.strobe = 1;

    lz77_divide:
        while (true) {
#pragma HLS LOOP_TRIPCOUNT min = 1048576 max = 1048576
#pragma HLS PIPELINE II = 1
            uint32_t ltreeIdx, dtreeIdx;
            if (next_state == WRITE_LITERAL) {
                auto encodedValue = inStream.read();
                if (encodedValue.strobe == 0) {
                    last_block = just_started;
                    just_started = true;
                    break;
                }
                just_started = false;
                tCh = encodedValue.data[0].range(7, 0);
                tLen = encodedValue.data[0].range(15, 8);
                uint16_t tOffset = encodedValue.data[0].range(31, 16);
                tOffset1 = encodedValue.data[0].range(23, 16);
                tOffset2 = encodedValue.data[0].range(31, 24);

                if (tLen) {
                    ltreeIdx = length_code[tLen] + 256 + 1;
                    dtreeIdx = d_code(tOffset, dist_code);
                    outValue.data[0].range(8, 8) = 1;
                    outValue.data[0].range(7, 0) = tLen;
                    next_state = WRITE_OFFSET0;
                } else {
                    ltreeIdx = tCh;
                    outValue.data[0] = tCh;
                    next_state = WRITE_LITERAL;
                }
            } else if (next_state == WRITE_OFFSET0) {
                outValue.data[0] = tOffset1;
                next_state = WRITE_OFFSET1;
            } else if (next_state == WRITE_OFFSET1) {
                outValue.data[0] = tOffset2;
                next_state = WRITE_LITERAL;
            }
            outStream << outValue;
        }
        // File/Block end Condition
        outValue.strobe = 0;
        outStream << outValue;
    }
}

/**
 * @brief Objective of this module is to pick character with
 * higher match length in the offset window range.
 *
 * @tparam MATCH_LEN match length
 * @tparam OFFSET_WINDOW offset window range
 *
 * @param inStream input stream
 * @param outStream output stream
 * @param input_size input stream size
 */
template <int MATCH_LEN, int OFFSET_WINDOW>
void lzBestMatchFilter(hls::stream<compressd_dt>& inStream, hls::stream<compressd_dt>& outStream, uint32_t input_size) {
    const int c_max_match_length = MATCH_LEN;
    if (input_size == 0) return;

    compressd_dt compare_window[MATCH_LEN];
#pragma HLS array_partition variable = compare_window

    // Initializing shift registers
    for (uint32_t i = 0; i < c_max_match_length; i++) {
#pragma HLS UNROLL
        compare_window[i] = inStream.read();
    }

lz_bestMatchFilter:
    for (uint32_t i = c_max_match_length; i < input_size; i++) {
#pragma HLS PIPELINE II = 1
        // shift register logic
        compressd_dt outValue = compare_window[0];
        for (uint32_t j = 0; j < c_max_match_length - 1; j++) {
#pragma HLS UNROLL
            compare_window[j] = compare_window[j + 1];
        }
        compare_window[c_max_match_length - 1] = inStream.read();

        uint8_t match_length = outValue.range(15, 8);
        bool best_match = 1;
        // Find Best match
        for (uint32_t j = 0; j < c_max_match_length; j++) {
            compressd_dt compareValue = compare_window[j];
            uint8_t compareLen = compareValue.range(15, 8);
            if (match_length + j < compareLen) {
                best_match = 0;
            }
        }
        if (best_match == 0) {
            outValue.range(15, 8) = 0;
            outValue.range(31, 16) = 0;
        }
        outStream << outValue;
    }

lz_bestMatchFilter_left_over:
    for (uint32_t i = 0; i < c_max_match_length; i++) {
#pragma HLS PIPELINE off
        outStream << compare_window[i];
    }
}

template <int MATCH_LEN, int OFFSET_WINDOW>
void lzBestMatchFilter(hls::stream<IntVectorStream_dt<32, 1> >& inStream,
                       hls::stream<IntVectorStream_dt<32, 1> >& outStream) {
    const int c_max_match_length = MATCH_LEN;

    IntVectorStream_dt<32, 1> compare_window;
    IntVectorStream_dt<32, 1> outStreamValue;

    while (true) {
        auto nextVal = inStream.read();
        if (nextVal.strobe == 0) {
            outStreamValue.strobe = 0;
            outStream << outStreamValue;
            break;
        }
        // assuming that, at least bytes more than LEFT_BYTES will be present at the input
        compare_window = nextVal;
        nextVal = inStream.read();

    lz_bestMatchFilter:
        for (; nextVal.strobe != 0;) {
#pragma HLS PIPELINE II = 1
            // shift register logic
            IntVectorStream_dt<32, 1> outValue = compare_window;
            compare_window = nextVal;
            nextVal = inStream.read();

            uint8_t match_length = outValue.data[0].range(15, 8);
            bool best_match = 1;
            // Find Best match
            compressd_dt compareValue = compare_window.data[0];
            uint8_t compareLen = compareValue.range(15, 8);
            if (match_length < compareLen) {
                best_match = 0;
            }
            if (best_match == 0) {
                outValue.data[0].range(15, 8) = 0;
                outValue.data[0].range(31, 16) = 0;
            }
            outStream << outValue;
        }
        outStream << compare_window;
        outStreamValue.strobe = 0;
        outStream << outStreamValue;
    }
}

/**
 * @brief This module helps in improving the compression ratio.
 * Finds a better match length by performing more character matches
 * with supported max match, while maintaining an offset window.
 * Booster offset Window template argument (default value is 16K)
 * internally consume BRAM memory to implement history window.
 * Higher the booster value can give better compression ratio but
 * will consume more BRAM resources.
 *
 * @tparam MAX_MATCH_LEN maximum length allowed for character match
 * @tparam BOOSTER_OFFSET_WINDOW offset window to store/match the character
 * @tparam LEFT_BYTES last 64 left over bytes
 *
 * @param inStream input stream 32bit per read
 * @param outStream output stream 32bit per write
 *
 */
template <int MAX_MATCH_LEN, int BLOCKSIZE = 16384, int BOOSTER_OFFSET_WINDOW = 8 * 1024, int LEFT_BYTES = 64>
void lzBooster(hls::stream<IntVectorStream_dt<32, 1> >& inStream, hls::stream<IntVectorStream_dt<32, 1> >& outStream) {
    constexpr uint16_t c_fifo_depth = LEFT_BYTES + 2;
    constexpr int c_boosterOffsetWindow = (BLOCKSIZE < BOOSTER_OFFSET_WINDOW) ? BLOCKSIZE : BOOSTER_OFFSET_WINDOW;
    constexpr int c_parallel_units = 2; // Reduced parallel units to save resources
    
    bool last_block = false;
    bool block_end = true;
    IntVectorStream_dt<32, 1> outStreamValue;
    ap_uint<32> outValue = 0;

    while (true) {
        // Single declaration of local memory - optimized size
        uint8_t local_mem[c_boosterOffsetWindow];
#pragma HLS ARRAY_PARTITION variable=local_mem cyclic factor=c_parallel_units dim=1
#pragma HLS BIND_STORAGE variable=local_mem type=RAM_1P impl=BRAM
        hls::stream<ap_uint<32> > lclBufStream("lclBufStream");
#pragma HLS STREAM variable = lclBufStream depth = c_fifo_depth
#pragma HLS BIND_STORAGE variable = lclBufStream type = fifo impl = srl
        uint32_t iIdx = 0;

        uint32_t match_loc = 0;
        uint32_t match_len = 0;

        outStreamValue.strobe = 1;

        bool matchFlag = false;
        bool outFlag = false;
        bool boostFlag = false;
        uint16_t skip_len = 0;
        
        // Pre-fetch match characters
        uint8_t nextMatchCh[c_parallel_units];
#pragma HLS ARRAY_PARTITION variable=nextMatchCh complete dim=1
        
        // Initialize match characters
        for (int i = 0; i < c_parallel_units; i++) {
#pragma HLS UNROLL
            nextMatchCh[i] = local_mem[(match_loc + i) % c_boosterOffsetWindow];
        }

        // check and exit if end of data
        auto nextVal = inStream.read();
        if (nextVal.strobe == 0) {
            outStreamValue.strobe = 0;
            outStream << outStreamValue;
            break;
        }

    // Buffer initialization with improved parallelism
    lz_booster_init_buf:
        for (uint16_t i = 0; i < LEFT_BYTES; i += c_parallel_units) {
#pragma HLS PIPELINE II = 1
            ap_uint<32> tmp[c_parallel_units];
#pragma HLS ARRAY_PARTITION variable=tmp complete dim=1
            
            // Read multiple values in parallel
            for (int j = 0; j < c_parallel_units && i + j < LEFT_BYTES; j++) {
#pragma HLS UNROLL
                if (j == 0) {
                    tmp[j] = nextVal.data[0];
                    nextVal = inStream.read();
                } else {
                    tmp[j] = nextVal.data[0];
                    nextVal = inStream.read();
                }
            }
            
            // Write to buffer in parallel
            for (int j = 0; j < c_parallel_units && i + j < LEFT_BYTES; j++) {
#pragma HLS UNROLL
                lclBufStream << tmp[j];
            }
        }

    lz_booster:
        for (; nextVal.strobe != 0; ++iIdx) { // iterate through block data
#pragma HLS PIPELINE II = 1
#pragma HLS LOOP_TRIPCOUNT min=1024 max=8192
#ifndef DISABLE_DEPENDENCE
#pragma HLS dependence variable = local_mem inter false
#pragma HLS dependence variable = nextMatchCh inter false
#endif
            // read value from fifo
            auto inValue = lclBufStream.read();
            // read from input stream
            auto tmp = inStream.read();
            // write new value to fifo
            lclBufStream << nextVal.data[0];
            // update next val
            nextVal = tmp;

            uint8_t tCh = inValue.range(7, 0);
            uint8_t tLen = inValue.range(15, 8);
            uint16_t tOffset = inValue.range(31, 16);

            if (tOffset < c_boosterOffsetWindow) {
                boostFlag = true;
            } else {
                boostFlag = false;
            }
            uint8_t match_ch = nextMatchCh[0]; // Use first element from the array
            local_mem[iIdx % c_boosterOffsetWindow] = tCh;
            outFlag = false;

            if (skip_len) {
                skip_len--;
            } else if (matchFlag && (match_len < MAX_MATCH_LEN) && (tCh == match_ch)) {
                match_len++;
                match_loc++;
                outValue.range(15, 8) = match_len;
            } else {
                match_len = 1;
                match_loc = iIdx - tOffset;
                if (iIdx) outFlag = true;
                outStreamValue.data[0] = outValue;
                outValue = inValue;
                if (tLen) {
                    if (boostFlag) {
                        matchFlag = true;
                        skip_len = 0;
                    } else {
                        matchFlag = false;
                        skip_len = tLen - 1;
                    }
                } else {
                    matchFlag = false;
                }
            }
            // Update all entries in nextMatchCh array
            for (int k = 0; k < c_parallel_units; k++) {
#pragma HLS UNROLL
                nextMatchCh[k] = local_mem[(match_loc + k) % c_boosterOffsetWindow];
            }
            if (outFlag) {
                outStream << outStreamValue;
            }
        }
        outStreamValue.data[0] = outValue;
        outStream << outStreamValue;
    // Balanced throughput processing for optimal performance
    constexpr unsigned int PACKING_FACTOR = 8;  // Balanced back to 8 for resource efficiency
    constexpr unsigned int FULL_PACKS = LEFT_BYTES / PACKING_FACTOR;
    constexpr unsigned int REMAINING_BYTES = LEFT_BYTES % PACKING_FACTOR;

    // Use balanced register arrays for good parallelism
    ap_uint<8> prefetch_buffer[PACKING_FACTOR];
#pragma HLS ARRAY_PARTITION variable=prefetch_buffer complete dim=1

    // Efficient prefetch - load 8 bytes in parallel
    for (uint8_t p = 0; p < PACKING_FACTOR; ++p) {
#pragma HLS UNROLL
        prefetch_buffer[p] = lclBufStream.read();
    }

lz_booster_left_bytes:
    for (uint16_t i = 0; i < FULL_PACKS - 1; ++i) {
#pragma HLS PIPELINE II = 1
#pragma HLS LOOP_TRIPCOUNT min=7 max=7  // LEFT_BYTES=64, PACKING_FACTOR=8

        // Use balanced data packing - process 2 32-bit words (8 bytes) per cycle
        ap_uint<32> packed_data1 = 0, packed_data2 = 0;
        
        // Balanced packing - 8 bytes in 2 words
        for (uint8_t p = 0; p < 4; ++p) {
#pragma HLS UNROLL
            packed_data1.range((p + 1) * 8 - 1, p * 8) = prefetch_buffer[p];
            packed_data2.range((p + 1) * 8 - 1, p * 8) = prefetch_buffer[p + 4];
        }

        // Pipeline-optimized prefetch for next iteration
        for (uint8_t p = 0; p < PACKING_FACTOR; ++p) {
#pragma HLS UNROLL
            prefetch_buffer[p] = lclBufStream.read();
        }

        // Parallel output of 2 32-bit words
        outStreamValue.data[0] = packed_data1;
        outStream << outStreamValue;
        outStreamValue.data[0] = packed_data2;
        outStream << outStreamValue;
    }

    // Process final complete pack efficiently
    ap_uint<32> final_packed1 = 0, final_packed2 = 0;
    for (uint8_t p = 0; p < 4; ++p) {
#pragma HLS UNROLL
        final_packed1.range((p + 1) * 8 - 1, p * 8) = prefetch_buffer[p];
        final_packed2.range((p + 1) * 8 - 1, p * 8) = prefetch_buffer[p + 4];
    }
    
    outStreamValue.data[0] = final_packed1;
    outStream << outStreamValue;
    outStreamValue.data[0] = final_packed2;
    outStream << outStreamValue;

    // 处理剩余字节（如果有）
    if (REMAINING_BYTES > 0) {
        ap_uint<32> remain_packed = 0;
#pragma HLS LATENCY min=1 max=1
        for (uint8_t p = 0; p < REMAINING_BYTES; ++p) {
#pragma HLS UNROLL
            ap_uint<8> current_byte = lclBufStream.read();
            remain_packed.range((p + 1) * 8 - 1, p * 8) = current_byte;
        }
        outStreamValue.data[0] = remain_packed;
        outStream << outStreamValue;
    }

    // Block end
        outStreamValue.strobe = 0;
        outStream << outStreamValue;
    }
}

/**
 * @brief This module helps in improving the compression ratio.
 * Finds a better match length by performing more character matches
 * with supported max match, while maintaining an offset window.
 * Booster offset Window template argument (default value is 16K)
 * internally consume BRAM memory to implement history window.
 * Higher the booster value can give better compression ratio but
 * will consume more BRAM resources.
 *
 * @tparam MAX_MATCH_LEN maximum length allowed for character match
 * @tparam BOOSTER_OFFSET_WINDOW offset window to store/match the character
 *
 * @param inStream input stream 32bit per read
 * @param outStream output stream 32bit per write
 * @param input_size input size
 * @param left_bytes last 64 left over bytes
 *
 */
template <int MAX_MATCH_LEN, int BOOSTER_OFFSET_WINDOW = 8 * 1024, int LEFT_BYTES = 64>
void lzBooster(hls::stream<compressd_dt>& inStream, hls::stream<compressd_dt>& outStream, uint32_t input_size) {
    if (input_size == 0) return;
    uint8_t local_mem[BOOSTER_OFFSET_WINDOW];
#pragma HLS BIND_STORAGE variable=local_mem type=RAM_1P impl=BRAM
    uint32_t match_loc = 0;
    uint32_t match_len = 0;
    compressd_dt outValue;
    compressd_dt outStreamValue;
    bool matchFlag = false;
    bool outFlag = false;
    bool boostFlag = false;
    uint16_t skip_len = 0;
    uint8_t nextMatchCh = local_mem[match_loc % BOOSTER_OFFSET_WINDOW];

lz_booster:
    for (uint32_t i = 0; i < (input_size - LEFT_BYTES); i++) {
#pragma HLS PIPELINE II = 1
#pragma HLS dependence variable = local_mem inter false
        compressd_dt inValue = inStream.read();
        uint8_t tCh = inValue.range(7, 0);
        uint8_t tLen = inValue.range(15, 8);
        uint16_t tOffset = inValue.range(31, 16);
        if (tOffset < BOOSTER_OFFSET_WINDOW) {
            boostFlag = true;
        } else {
            boostFlag = false;
        }
        uint8_t match_ch = nextMatchCh;
        local_mem[i % BOOSTER_OFFSET_WINDOW] = tCh;
        outFlag = false;

        if (skip_len) {
            skip_len--;
        } else if (matchFlag && (match_len < MAX_MATCH_LEN) && (tCh == match_ch)) {
            match_len++;
            match_loc++;
            outValue.range(15, 8) = match_len;
        } else {
            match_len = 1;
            match_loc = i - tOffset;
            if (i) outFlag = true;
            outStreamValue = outValue;
            outValue = inValue;
            if (tLen) {
                if (boostFlag) {
                    matchFlag = true;
                    skip_len = 0;
                } else {
                    matchFlag = false;
                    skip_len = tLen - 1;
                }
            } else {
                matchFlag = false;
            }
        }
        nextMatchCh = local_mem[match_loc % BOOSTER_OFFSET_WINDOW];
        if (outFlag) outStream << outStreamValue;
    }
    outStream << outValue;

    // 定义高效的并行处理参数
    constexpr unsigned int PACKING_FACTOR = 16;  // 每次处理16字节
    constexpr unsigned int WORDS_PER_PACK = PACKING_FACTOR / 4;  // 每组4个32位字
    constexpr unsigned int FULL_PACKS = LEFT_BYTES / PACKING_FACTOR;

    // 优化的缓存结构
    struct {
        ap_uint<8> bytes[PACKING_FACTOR];
        ap_uint<32> words[WORDS_PER_PACK];
    } buffer;
#pragma HLS ARRAY_PARTITION variable=buffer.bytes complete dim=1
#pragma HLS ARRAY_PARTITION variable=buffer.words complete dim=1

lz_booster_left_bytes:
    for (uint16_t i = 0; i < FULL_PACKS; ++i) {
#pragma HLS PIPELINE II = 1
#pragma HLS LOOP_TRIPCOUNT min=4 max=4
        
        // 并行读取和处理4字节组
        for (uint8_t w = 0; w < WORDS_PER_PACK; ++w) {
#pragma HLS UNROLL
            compressd_dt packed_word = 0;
            
            // 读取并打包4个字节
            for (uint8_t b = 0; b < 4; ++b) {
#pragma HLS UNROLL
                compressd_dt temp = inStream.read();
                packed_word.range((b + 1) * 8 - 1, b * 8) = temp.range(7, 0);
            }
            
            outStream << packed_word;
        }
    }

    // 处理剩余字节
    for (uint16_t i = FULL_PACKS * PACKING_FACTOR; i < LEFT_BYTES; ++i) {
#pragma HLS PIPELINE II = 1
        outStreamValue = inStream.read();
        outStream << outStreamValue;
    }
}

/**
 * @brief This module checks if match length exists, and if
 * match length exists it filters the match length -1 characters
 * writing to output stream.
 *
 * @tparam MATCH_LEN length of matched segment
 * @tparam OFFSET_WINDOW output window
 *
 * @param inStream input stream
 * @param outStream output stream
 * @param input_size input stream size
 * @param left_bytes bytes left in block
 */
template <int LEFT_BYTES = 64>
static void lzFilter(hls::stream<compressd_dt>& inStream, hls::stream<compressd_dt>& outStream, uint32_t input_size) {
    if (input_size == 0) return;
    uint32_t skip_len = 0;
lz_filter:
    for (uint32_t i = 0; i < input_size - LEFT_BYTES; i++) {
#pragma HLS PIPELINE II = 1
        compressd_dt inValue = inStream.read();
        uint8_t tLen = inValue.range(15, 8);
        if (skip_len) {
            skip_len--;
        } else {
            outStream << inValue;
            if (tLen) skip_len = tLen - 1;
        }
    }
lz_filter_left_bytes:
    for (uint32_t i = 0; i < LEFT_BYTES; i++) {
        outStream << inStream.read();
    }
}

} // namespace compression
} // namespace xf
#endif // _XFCOMPRESSION_LZ_OPTIONAL_HPP_
