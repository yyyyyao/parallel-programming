/*
 * =====================================================================================
 *
 *       Filename:  blur.cpp
 *
 *     
 *
 *        Version:  1.0
 *        Created:  2013/07/28 14時20分48秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tetsuhiko Yao (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <iostream>
#include <string>
#include <stdio.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#define int8 char
#define uint8 unsigned char
#define uint16 unsigned short
#define uint32 unsigned int
#define uint unsigned int

// Various JPEG enums and tables.
enum { M_SOF0 = 0xC0, M_DHT = 0xC4, M_SOI = 0xD8, M_EOI = 0xD9, M_SOS = 0xDA, M_DQT = 0xDB, M_APP0 = 0xE0 };
enum { DC_LUM_CODES = 12, AC_LUM_CODES = 256, DC_CHROMA_CODES = 12, AC_CHROMA_CODES = 256, MAX_HUFF_SYMBOLS = 257, MAX_HUFF_CODESIZE = 32 };

static uint8 s_dc_lum_bits[17] = { 0,0,1,5,1,1,1,1,1,1,0,0,0,0,0,0,0 };
static uint8 s_dc_lum_val[DC_LUM_CODES] = { 0,1,2,3,4,5,6,7,8,9,10,11 };
static uint8 s_ac_lum_bits[17] = { 0,0,2,1,3,3,2,4,3,5,5,4,4,0,0,1,0x7d };
static uint8 s_ac_lum_val[AC_LUM_CODES]  =
{
  0x01,0x02,0x03,0x00,0x04,0x11,0x05,0x12,0x21,0x31,0x41,0x06,0x13,0x51,0x61,0x07,0x22,0x71,0x14,0x32,0x81,0x91,0xa1,0x08,0x23,0x42,0xb1,0xc1,0x15,0x52,0xd1,0xf0,
  0x24,0x33,0x62,0x72,0x82,0x09,0x0a,0x16,0x17,0x18,0x19,0x1a,0x25,0x26,0x27,0x28,0x29,0x2a,0x34,0x35,0x36,0x37,0x38,0x39,0x3a,0x43,0x44,0x45,0x46,0x47,0x48,0x49,
  0x4a,0x53,0x54,0x55,0x56,0x57,0x58,0x59,0x5a,0x63,0x64,0x65,0x66,0x67,0x68,0x69,0x6a,0x73,0x74,0x75,0x76,0x77,0x78,0x79,0x7a,0x83,0x84,0x85,0x86,0x87,0x88,0x89,
  0x8a,0x92,0x93,0x94,0x95,0x96,0x97,0x98,0x99,0x9a,0xa2,0xa3,0xa4,0xa5,0xa6,0xa7,0xa8,0xa9,0xaa,0xb2,0xb3,0xb4,0xb5,0xb6,0xb7,0xb8,0xb9,0xba,0xc2,0xc3,0xc4,0xc5,
  0xc6,0xc7,0xc8,0xc9,0xca,0xd2,0xd3,0xd4,0xd5,0xd6,0xd7,0xd8,0xd9,0xda,0xe1,0xe2,0xe3,0xe4,0xe5,0xe6,0xe7,0xe8,0xe9,0xea,0xf1,0xf2,0xf3,0xf4,0xf5,0xf6,0xf7,0xf8,
  0xf9,0xfa
};
static uint8 s_dc_chroma_bits[17] = { 0,0,3,1,1,1,1,1,1,1,1,1,0,0,0,0,0 };
static uint8 s_dc_chroma_val[DC_CHROMA_CODES]  = { 0,1,2,3,4,5,6,7,8,9,10,11 };
static uint8 s_ac_chroma_bits[17] = { 0,0,2,1,2,4,4,3,4,7,5,4,4,0,1,2,0x77 };
static uint8 s_ac_chroma_val[AC_CHROMA_CODES] =
{
  0x00,0x01,0x02,0x03,0x11,0x04,0x05,0x21,0x31,0x06,0x12,0x41,0x51,0x07,0x61,0x71,0x13,0x22,0x32,0x81,0x08,0x14,0x42,0x91,0xa1,0xb1,0xc1,0x09,0x23,0x33,0x52,0xf0,
  0x15,0x62,0x72,0xd1,0x0a,0x16,0x24,0x34,0xe1,0x25,0xf1,0x17,0x18,0x19,0x1a,0x26,0x27,0x28,0x29,0x2a,0x35,0x36,0x37,0x38,0x39,0x3a,0x43,0x44,0x45,0x46,0x47,0x48,
  0x49,0x4a,0x53,0x54,0x55,0x56,0x57,0x58,0x59,0x5a,0x63,0x64,0x65,0x66,0x67,0x68,0x69,0x6a,0x73,0x74,0x75,0x76,0x77,0x78,0x79,0x7a,0x82,0x83,0x84,0x85,0x86,0x87,
  0x88,0x89,0x8a,0x92,0x93,0x94,0x95,0x96,0x97,0x98,0x99,0x9a,0xa2,0xa3,0xa4,0xa5,0xa6,0xa7,0xa8,0xa9,0xaa,0xb2,0xb3,0xb4,0xb5,0xb6,0xb7,0xb8,0xb9,0xba,0xc2,0xc3,
  0xc4,0xc5,0xc6,0xc7,0xc8,0xc9,0xca,0xd2,0xd3,0xd4,0xd5,0xd6,0xd7,0xd8,0xd9,0xda,0xe2,0xe3,0xe4,0xe5,0xe6,0xe7,0xe8,0xe9,0xea,0xf2,0xf3,0xf4,0xf5,0xf6,0xf7,0xf8,
  0xf9,0xfa
};

static uint8 s_zag[64] = { 0,1,8,16,9,2,3,10,17,24,32,25,18,11,4,5,12,19,26,33,40,48,41,34,27,20,13,6,7,14,21,28,35,42,49,56,57,50,43,36,29,22,15,23,30,37,44,51,58,59,52,45,38,31,39,46,53,60,61,54,47,55,62,63 };

static int s_std_lum_quant[64] = { 16,11,12,14,12,10,16,14,13,14,18,17,16,19,24,40,26,24,22,22,24,49,35,37,29,40,58,51,61,60,57,51,56,55,64,72,92,78,64,68,87,69,55,56,80,109,81,87,95,98,103,104,103,62,77,113,121,112,100,120,92,101,103,99 };
static int s_std_croma_quant[64] = { 17,18,18,24,21,24,47,26,26,47,99,66,56,66,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99 };

const int YR = 19595, YG = 38470, YB = 7471, CB_R = -11059, CB_G = -21709, CB_B = 32768, CR_R = 32768, CR_G = -27439, CR_B = -5329;
static inline uint8 clamp(int i) { if (static_cast<uint>(i) > 255U) { if (i < 0) i = 0; else if (i > 255) i = 255; } return static_cast<uint8>(i); }

/* variables.these variables should be c++ member variables. */
int m_bits_in;
int m_num_components = 3;
int* m_quantization_tables[3];
int m_image_y, m_image_x;

uint8 m_comp_h_samp[3] = {2, 1, 1};
uint8 m_comp_v_samp[3] = {2, 1, 1};
/*
m_comp_h_samp[0] = 2; m_comp_v_samp[0] = 2;
m_comp_h_samp[1] = 1; m_comp_v_samp[1] = 1;
m_comp_h_samp[2] = 1; m_comp_v_samp[2] = 1;
*/

#define WRITE_FILES 1
#ifdef WRITE_FILES
/* for writing bits variables */
#define JPGE_OUT_BUF_SIZE 2048
unsigned int m_out_buf_left = JPGE_OUT_BUF_SIZE;
uint8 m_out_buf[JPGE_OUT_BUF_SIZE];
uint8* m_pOut_buf = m_out_buf;
bool m_all_stream_writes_succeeded;

unsigned int m_bit_buffer = 0;
unsigned int m_bit_in;
bool m_bStatus = true;
FILE* m_pFile;
/* for writing bits functions */

bool close()
{
  if (m_pFile)
  {
     if (fclose(m_pFile) == EOF)
     {
        m_bStatus = false;
     }
     m_pFile = NULL;
  }
  return m_bStatus;
}

bool open(const char *pFilename)
{
      close();
      m_pFile = fopen(pFilename, "wb");
      m_bStatus = (m_pFile != NULL);
      return m_bStatus;
}

bool put_buf(const void* pBuf, int len)
{
  m_bStatus = m_bStatus && (fwrite(pBuf, len, 1, m_pFile) == 1);
  return m_bStatus;
}

template<class T> inline bool put_obj(const T& obj) { return put_buf(&obj, sizeof(T)); }


// JPEG marker generation.
#if 0
void emit_byte(uint8 i)
{
  m_all_stream_writes_succeeded = m_all_stream_writes_succeeded && put_obj(i);
}
#else
void emit_byte(uint8 i)
{
  bool temp = put_obj(i);
  m_all_stream_writes_succeeded = m_all_stream_writes_succeeded && temp;
}
#endif

void emit_word(uint i)
{
  emit_byte((uint8)(i >> 8)); emit_byte((uint8)(i & 0xFF));
}

void emit_marker(int marker)
{
  emit_byte((int8)(0xFF)); emit_byte((uint8)(marker));
}

// Emit JFIF marker
void emit_jfif_app0()
{
  emit_marker(M_APP0);
  emit_word(2 + 4 + 1 + 2 + 1 + 2 + 2 + 1 + 1);
  emit_byte(0x4A); emit_byte(0x46); emit_byte(0x49); emit_byte(0x46); /* Identifier: ASCII "JFIF" */
  emit_byte(0);
  emit_byte(1);      /* Major version */
  emit_byte(1);      /* Minor version */
  emit_byte(0);      /* Density unit */
  emit_word(1);
  emit_word(1);
  emit_byte(0);      /* No thumbnail image */
  emit_byte(0);
}

// Emit quantization tables
void emit_dqt()
{
  for (int i = 0; i < ((m_num_components == 3) ? 2 : 1); i++)
  {
    emit_marker(M_DQT);
    emit_word(64 + 1 + 2);
    emit_byte(static_cast<uint8>(i));
    for (int j = 0; j < 64; j++)
      emit_byte(static_cast<uint8>(m_quantization_tables[i][j]));
  }
}

// Emit start of frame marker
void emit_sof()
{
  emit_marker(M_SOF0);                           /* baseline */
  emit_word(3 * m_num_components + 2 + 5 + 1);
  emit_byte(8);                                  /* precision */
  emit_word(m_image_y);
  emit_word(m_image_x);
  emit_byte(m_num_components);
  for (int i = 0; i < m_num_components; i++)
  {
    emit_byte(static_cast<uint8>(i + 1));                                   /* component ID     */
    emit_byte((m_comp_h_samp[i] << 4) + m_comp_v_samp[i]);  /* h and v sampling */
    emit_byte(i > 0);                                   /* quant. table num */
  }
}

// Emit Huffman table.
void emit_dht(uint8 *bits, uint8 *val, int index, bool ac_flag)
{
  emit_marker(M_DHT);

  int length = 0;
  for (int i = 1; i <= 16; i++)
    length += bits[i];

  emit_word(length + 2 + 1 + 16);
  emit_byte(static_cast<uint8>(index + (ac_flag << 4)));

  for (int i = 1; i <= 16; i++)
    emit_byte(bits[i]);

  for (int i = 0; i < length; i++)
    emit_byte(val[i]);
}

// Emit all Huffman tables.
void emit_dhts()
{
  emit_dht(s_dc_lum_bits, s_dc_lum_val, 0, false);
  emit_dht(s_ac_lum_bits, s_ac_lum_val, 0, true);
  if (m_num_components == 3)
  {
    emit_dht(s_dc_chroma_bits, s_dc_chroma_val, 1, false);
    emit_dht(s_ac_chroma_bits, s_ac_chroma_val, 1, true);
  }
}

// emit start of scan
void emit_sos()
{
  emit_marker(M_SOS);
  emit_word(2 * m_num_components + 2 + 1 + 3);
  emit_byte(m_num_components);
  for (int i = 0; i < m_num_components; i++)
  {
    emit_byte(static_cast<uint8>(i + 1));
    if (i == 0)
      emit_byte((0 << 4) + 0);
    else
      emit_byte((1 << 4) + 1);
  }
  emit_byte(0);     /* spectral selection */
  emit_byte(63);
  emit_byte(0);
}

// Emit all markers at beginning of image file.
void emit_markers()
{
  emit_marker(M_SOI);
  emit_jfif_app0();
  emit_dqt();
  emit_sof();
  emit_dhts();
  emit_sos();

}

#if 1
void flush_output_buffer()
{
  if (m_out_buf_left != JPGE_OUT_BUF_SIZE) {
    bool temp;
    temp = put_buf(m_out_buf, JPGE_OUT_BUF_SIZE - m_out_buf_left);
    m_all_stream_writes_succeeded = m_all_stream_writes_succeeded && temp;
  }
  m_pOut_buf = m_out_buf;
  m_out_buf_left = JPGE_OUT_BUF_SIZE;
}
#else
void flush_output_buffer()
{
  if (m_out_buf_left != JPGE_OUT_BUF_SIZE)
    m_all_stream_writes_succeeded = m_all_stream_writes_succeeded && put_buf(m_out_buf, JPGE_OUT_BUF_SIZE - m_out_buf_left);
  m_pOut_buf = m_out_buf;
  m_out_buf_left = JPGE_OUT_BUF_SIZE;
}
#endif

#if 1
#define JPGE_PUT_BYTE(c) { *m_pOut_buf++ = (c); if (--m_out_buf_left == 0) flush_output_buffer(); } /* m_out_buf_left is 2048. */

#else
void JPGE_PUT_BYTE(uint8 c) {
  *m_pOut_buf++ = (c);
  if (--m_out_buf_left == 0)
    flush_output_buffer();
}
#endif

void put_bits(uint bits, uint len) /* uint is unsigned int (=32 bit) */
{
  m_bit_buffer |= ((uint32)bits << (24 - (m_bits_in += len))); /* m_bits_in shows the data length of buffer. m_bits_in calc immediately.*/
  while (m_bits_in >= 8)
  {
    uint8 c;

    JPGE_PUT_BYTE(c = (uint8)((m_bit_buffer >> 16) & 0xFF));
    if (c == 0xFF) JPGE_PUT_BYTE(0); /* if bit is 0XFF, write zero? why? */
    m_bit_buffer <<= 8;
    m_bits_in -= 8;
  }
}
bool terminate_pass_two()
{
  //put_bits(0x7F, 7);
  put_bits(0x7F, 8);
  flush_output_buffer();
  emit_marker(M_EOI);
  return true;
}
#endif /* WRITE_FILES */

void compute_huffman_table(uint *codes, uint8 *code_sizes, const uint8 *bits, const uint8 *val)
{
  int i, l, last_p, si;
  uint8 huff_size[257];
  uint huff_code[257];
  uint code;

  int p = 0;
  for (l = 1; l <= 16; l++) /* 2^16 = 65536? */
    for (i = 1; i <= bits[l]; i++)
      huff_size[p++] = (char)l; /* huff_size is the length of huffman mark */

  huff_size[p] = 0; last_p = p; // write sentinel      stopper?

  code = 0; si = huff_size[0]; p = 0;

  while (huff_size[p]) /* until reach the stopper.*/
  {
    while (huff_size[p] == si)
      huff_code[p++] = code++;
    code <<= 1;
    si++;
  }

  memset(codes, 0, sizeof(codes[0])*256);
  memset(code_sizes, 0, sizeof(code_sizes[0])*256);
  for (p = 0; p < last_p; p++)
  {
    codes[val[p]]      = huff_code[p];
    code_sizes[val[p]] = huff_size[p];
    printf("size:%d\tcode:%d\n", huff_size[p], huff_code[p]);
  }
  printf("-------------------------------\n");
}

void code_coefficients_pass_two(int* m_coefficient_array, const uint* dc_codes, const uint8* dc_code_sizes,
    const uint* ac_codes, const uint8* ac_code_sizes, int* m_last_dc_val )
{
  int i, j, run_len, nbits, temp1, temp2;
  int *pSrc = m_coefficient_array; /* quantized values. */
  const uint *codes[2]; /*  codes[0] is huffman code of dc value. codes[1] is av value. */
  const uint8 *code_sizes[2];

  codes[0] = dc_codes;
  code_sizes[0] = dc_code_sizes;
  codes[1] = ac_codes;
  code_sizes[1] = ac_code_sizes;
#if 0
  if (component_num == 0) /* equals to Y. */
  {
    codes[0] = m_huff_codes[0 + 0]; codes[1] = m_huff_codes[2 + 0]; /* codes[0] has address of lum dc huffman codes. codes[1] also lum ac.*/
    code_sizes[0] = m_huff_code_sizes[0 + 0]; code_sizes[1] = m_huff_code_sizes[2 + 0]; /* also has the huffman code sizes. */
  }
  else /* Cb or Cr. */
  {
    codes[0] = m_huff_codes[0 + 1]; codes[1] = m_huff_codes[2 + 1];
    code_sizes[0] = m_huff_code_sizes[0 + 1]; code_sizes[1] = m_huff_code_sizes[2 + 1];
  }
#endif

  temp1 = temp2 = pSrc[0] - *m_last_dc_val; /* minus by last dc value. */
  *m_last_dc_val = pSrc[0]; /* change to new last dc value. */

  if (temp1 < 0)
  {
    temp1 = -temp1; temp2--; /* temp1 is the diff former dc value and present dc value. */
  }

  nbits = 0;
  while (temp1)
  {
    nbits++; temp1 >>= 1; /* divided by 2 nbits shows bit length of diffs */
  }

  put_bits(codes[0][nbits], code_sizes[0][nbits]); /* puts the huffman codes. this codes shows data-length. */
  if (nbits) put_bits(temp2 & ((1 << nbits) - 1), nbits); /* if data is zero, we don't have to write data. in other case, we write the data.*/

  for (run_len = 0, i = 1; i < 64; i++) /* this loop handles ac values. */
  {
    if ((temp1 = m_coefficient_array[i]) == 0) /* if quantized value is 0 */
      run_len++; /* inc run-length. */
    else
    {
      while (run_len >= 16) /* if run-length exceeds 16, we have to put ZRL marks.*/
      {
        put_bits(codes[1][0xF0], code_sizes[1][0xF0]); /* put ZRL. */
        run_len -= 16;
      }
      if ((temp2 = temp1) < 0)
      {
        temp1 = -temp1;
        temp2--;
      }
      nbits = 1;
      while (temp1 >>= 1) /* if not zero data comes. */
        nbits++;
      j = (run_len << 4) + nbits; /* ac huffman code shows 0-run length and length of data size. */
      put_bits(codes[1][j], code_sizes[1][j]); /* put huffman codes. */
      put_bits(temp2 & ((1 << nbits) - 1), nbits); /* put data */
      run_len = 0;
    }
  }
  if (run_len) /* if run-length is remain, but reach the index of(8, 8) we should put EOB */
    put_bits(codes[1][0], code_sizes[1][0]);
}

void zigZagScan(const int* pSrc, int* pDst)
{
  for(int i = 0; i < 64; i++) {
    pDst[i] = pSrc[s_zag[i]];
  }
}
void applyQuantTable(const int* pSrc, int* pDst, int* quantTable)
{
  for(int i = 0; i < 8; i++) {
    for(int j = 0; j < 8; j++) {
      pDst[i * 8 + j] = pSrc[i * 8 + j] / quantTable[i * 8 + j];
    }
  }
}

void calcQuantTableByQuality(int* pSrc, int* pDst, int quality)
{
    int s = (quality < 50) ? (5000 / quality) : (200 - (2 * quality));
    for ( int i = 0; i < 64; i++ ) {
        int value = (s * pSrc[i] + 50) / 100;
        if ( value == 0 ) {
            value = 1;
        }
        if ( value > 255 ) {
            value = 255;
        }
        pDst[i] = value;
    }
}
void centerization(const uint8* pSrc, int8* pDst, int width, int height)
{
  for(int i = 0; i < height; i++) {
    for(int j = 0; j < width; j++) {
      pDst[i * width + j] = (int8)(pSrc[i * width + j] - 128);
    }
  }

}
void twoDimDct(const int8* pSrc, int* pDst)
{
  double cu, cv;
  double xCof, yCof;
  double cucv, val;
  for(int v = 0;v < 8; v++) {
    if(v) cv = 1;
    else cv = 1 / sqrt(2);
    for(int u = 0; u < 8; u++) {
      val = 0;
      if(u) cu = 1;
      else cu = 1 / sqrt(2);
      for(int y = 0; y < 8; y++) {
        yCof = (2 * y + 1) * v * M_PI / (2 * 8);
        for(int x = 0; x < 8; x++) {
          xCof = (2 * x + 1) * u * M_PI / (2 * 8);
          val += (int)(pSrc[y * 8 + x]) * cos(xCof) * cos(yCof);
        }
      }
      cucv = (2 * cu * cv) / 8;
      pDst[v * 8 + u] = floor(val * cucv + 0.5);
    }
  }
}
void reverseBlockSplitting(const uint8* pSrc, uint8* pDst, int width, int height)
{
  const int blockWidth = width / 8;
  for(int y = 0; y < height; y++) {
    for(int x = 0; x < width; x++) {
      int dstIndex = ((x % width) / 8 + (y / 8) * blockWidth) * 64 +
        (y % 8) * 8 + x % 8;
      int srcIndex = x + y * width;
      pDst[srcIndex] = pSrc[dstIndex];
    }
  }
}
void blockSplitMcuLine(const uint8* pSrc, uint8* pDst, int width)
{
  for(int i = 0; i < 8; i++) { /* because MUC height is 8. */
    for(int j = 0; j < width; j++) {
      int insideMcuIndex = i * 8 + j % 8;
      int mcuNum = j / 8;
      pDst[mcuNum * 64 + insideMcuIndex] = pSrc[i * width + j];
    }
  }
}
void SubSamplingH2V2(const uint8* pSrc, uint8* pDst, int width, int height)
{
  int halfHeight = ceil(height / 2);
  int halfWidth = ceil(width / 2);

  for(int i = 0; i < halfHeight; i++) {
    for(int j = 0; j < halfWidth; j++) {
      pDst[i * halfWidth + j] = pSrc[(i * 2) * width + (j * 2)];
    }
  }
}
void deSplitFromComp(uint8* pSrc, const uint8* pY, const uint8* pCb, const uint8* pCr,
    int height, int width)
{
  for(int i = 0; i < height; i++) {
    for(int j = 0; j < width; j++) {
      pSrc[i * width * 3 + j * 3 + 0] = pY[i * width + j] ;
      pSrc[i * width * 3 + j * 3 + 1] = pCb[i * width + j];
      pSrc[i * width * 3 + j * 3 + 2] = pCr[i * width + j];
    }
  }
}

void splitToComp(const uint8* pSrc, uint8* pY, uint8* pCb, uint8* pCr, int height, int width)
{
  for(int i = 0; i < height; i++) {
    for(int j = 0; j < width; j++) {
      pY[i * width + j] = pSrc[i * width * 3 + j * 3 + 0];
      pCb[i * width + j] = pSrc[i * width * 3 + j * 3 + 1];
      pCr[i * width + j] = pSrc[i * width * 3 + j * 3 + 2];
    }
  }
}

static void RGB_to_YCC(uint8* pDst, const uint8 *pSrc, int num_pixels)
{
  for ( ; num_pixels; pDst += 3, pSrc += 3, num_pixels--)
  {
    const int r = pSrc[2], g = pSrc[1], b = pSrc[0];
    pDst[0] = static_cast<uint8>((r * YR + g * YG + b * YB + 32768) >> 16);
    pDst[1] = clamp(128 + ((r * CB_R + g * CB_G + b * CB_B + 32768) >> 16));
    pDst[2] = clamp(128 + ((r * CR_R + g * CR_G + b * CR_B + 32768) >> 16));
  }
}

using namespace std;
int main(int argc, char* argv[])
{
  string windowName = "result window";
  cv::Mat input;
  cv::Mat refInput;
  cv::Size ksize;
 
  if(argc > 1) {
    input = cv::imread(argv[1], 1);
    refInput = cv::imread(argv[1], 0);
  }
  int height = input.rows;
  int width = input.cols;
  int halfWidth = ceil(width / 2);
  int halfHeight = ceil(height / 2);
  int comp = input.elemSize();

  printf("height:%d width:%d comp:%d\n", height, width, comp);

  cv::Mat output(input.rows, input.cols, CV_8UC3);
  cv::Mat grayOutput(input.rows, input.cols, CV_8UC1);

  ksize.width = 3;
  ksize.height = 3;
  if(argc > 2) {
    ksize.width = atoi(argv[2]);
    ksize.height = atoi(argv[2]);
  }

#ifdef WRITE_FILES
  open("output.jpg");
#endif
  unsigned char* srcPtr = (unsigned char*)input.data;
  unsigned char* refPtr = (unsigned char*)refInput.data;
  unsigned char* dstPtr = (unsigned char*)output.data;
  unsigned char* gDstPtr = (unsigned char*)grayOutput.data;
  unsigned char* srcPtrBuf;
  unsigned char* dstPtrBuf;

  /* at first, make tables. */
  /* Generates quantization tables */
  int quantTableY[64];
  int quantTableCbCr[64];
  int quality = 90;
  m_quantization_tables[0] = quantTableY;
  m_quantization_tables[1] = quantTableCbCr;
  m_quantization_tables[2] = quantTableCbCr;
  calcQuantTableByQuality(s_std_lum_quant, quantTableY, quality);
  calcQuantTableByQuality(s_std_croma_quant, quantTableCbCr, quality);

  /* Generates huffman codes */
  unsigned int huffCodesYDc[256];
  unsigned int huffCodesYAc[256];
  unsigned int huffCodesCbCrDc[256];
  unsigned int huffCodesCbCrAc[256];
  uint8 huffCodesSizesYDc[256];
  uint8 huffCodesSizesYAc[256];
  uint8 huffCodesSizesCbCrDc[256];
  uint8 huffCodesSizesCbCrAc[256];

  /* compute huffman Y Table. */
  compute_huffman_table(huffCodesYDc, huffCodesSizesYDc, s_dc_lum_bits, s_dc_lum_val);
  compute_huffman_table(huffCodesYAc, huffCodesSizesYAc, s_ac_lum_bits, s_ac_lum_val);
  /* compute huffman CbCr Table. */
  compute_huffman_table(huffCodesCbCrDc, huffCodesSizesCbCrDc, s_dc_chroma_bits, s_dc_chroma_val);
  compute_huffman_table(huffCodesCbCrAc, huffCodesSizesCbCrAc, s_ac_chroma_bits, s_ac_chroma_val);

  emit_markers();

  /* execute refrence codes */

  /*  RGB to YCbCr. */
  for(int i = 0; i < height; i++) {
    srcPtrBuf = srcPtr + i * width * comp;
    dstPtrBuf = dstPtr + i * width * comp;
    RGB_to_YCC(dstPtrBuf, srcPtrBuf, width);
  }

  /* split data for each comp. */
  cv::Mat Y(input.rows, input.cols, CV_8UC1);
  unsigned char* YPtr = (unsigned char*)Y.data;
  cv::Mat Cb(input.rows, input.cols, CV_8UC1);
  unsigned char* CbPtr = (unsigned char*)Cb.data;
  cv::Mat Cr(input.rows, input.cols, CV_8UC1);
  unsigned char* CrPtr = (unsigned char*)Cr.data;

  splitToComp(dstPtr, YPtr, CbPtr, CrPtr, width, height);

  /* SubSampling H2V2 */
  cv::Mat subCb(halfHeight, halfWidth, CV_8UC1);
  unsigned char* subCbPtr = (unsigned char*)subCb.data;
  cv::Mat subCr(halfHeight, halfWidth, CV_8UC1);
  unsigned char* subCrPtr = (unsigned char*)subCr.data;
  SubSamplingH2V2(CbPtr, subCbPtr, width, height);
  SubSamplingH2V2(CrPtr, subCrPtr, width, height);

  cv::Mat blockY(input.rows, input.cols, CV_8UC1);
  unsigned char* blockYPtr = (unsigned char*)blockY.data;
  unsigned char* YPtrBuf;
  unsigned char* blockYPtrBuf;
  /* block Splitting */
  for(int i = 0; i < height; i+=8) {
    YPtrBuf = YPtr + i * width ;
    blockYPtrBuf = blockYPtr + i * width;
    blockSplitMcuLine(YPtrBuf, blockYPtrBuf, width);
  }

  cv::Mat blockCb(halfHeight, halfWidth, CV_8UC1);
  unsigned char* blockCbPtr = (unsigned char*)blockCb.data;
  unsigned char* subCrPtrBuf;
  unsigned char* blockCrPtrBuf;
  cv::Mat blockCr(halfHeight, halfWidth, CV_8UC1);
  unsigned char* blockCrPtr = (unsigned char*)blockCr.data;
  unsigned char* subCbPtrBuf;
  unsigned char* blockCbPtrBuf;

  /* block Splitting */
  for(int i = 0; i < halfHeight; i+=8) {
    subCbPtrBuf = subCbPtr + i * halfWidth ;
    blockCbPtrBuf = blockCbPtr + i * halfWidth;
    blockSplitMcuLine(subCbPtrBuf, blockCbPtrBuf, halfWidth);

    subCrPtrBuf = subCrPtr + i * halfWidth ;
    blockCrPtrBuf = blockCrPtr + i * halfWidth;
    blockSplitMcuLine(subCrPtrBuf, blockCrPtrBuf, halfWidth);
  }

  /* apply reverse block splitting */
  cv::Mat reverseBlockY(input.rows, input.cols, CV_8UC1);
  unsigned char* reverseBlockYPtr = (unsigned char*)reverseBlockY.data;
  reverseBlockSplitting(blockYPtr, reverseBlockYPtr, width, height);

  cv::Mat reverseBlockCb(halfHeight, halfWidth, CV_8UC1);
  unsigned char* reverseBlockCbPtr = (unsigned char*)reverseBlockCb.data;
  reverseBlockSplitting(blockCbPtr, reverseBlockCbPtr, halfWidth, halfHeight);

  /* Centerization for DCT */
  cv::Mat centerY(input.rows, input.cols, CV_8UC1);
  char* centerYPtr= (char*)centerY.data;
  centerization(blockYPtr, centerYPtr, width, height);

  cv::Mat centerCb(halfHeight, halfWidth, CV_8UC1);
  char* centerCbPtr= (char*)centerCb.data;
  centerization(blockCbPtr, centerCbPtr, halfWidth, halfHeight);

  cv::Mat centerCr(halfHeight, halfWidth, CV_8UC1);
  char* centerCrPtr= (char*)centerCr.data;
  centerization(blockCrPtr, centerCrPtr, halfWidth, halfHeight);

  /* 2D DCT */
  cv::Mat dctY(input.rows, input.cols, CV_32SC1);
  int* dctYPtr= (int*)dctY.data;
  cv::Mat dctCr(halfHeight, halfWidth, CV_32SC1);
  int* dctCrPtr= (int*)dctCr.data;
  cv::Mat dctCb(halfHeight, halfWidth, CV_32SC1);
  int* dctCbPtr= (int*)dctCb.data;

  char* centerYPtrBuf;
  char* centerCbPtrBuf;
  char* centerCrPtrBuf;
  int* dctYPtrBuf;
  int* dctCbPtrBuf;
  int* dctCrPtrBuf;

  for(int i = 0; i < height; i+=8) {
    for(int j = 0; j < width; j+=8) {
      centerYPtrBuf = centerYPtr + i * width + j;
      dctYPtrBuf = dctYPtr + i * width + j;
      twoDimDct(centerYPtrBuf, dctYPtrBuf);
    }
  }

  for(int i = 0; i < halfHeight; i+=8) {
    for(int j = 0; j < halfWidth; j+=8) {
      centerCbPtrBuf = centerCbPtr + i * halfWidth + j;
      dctCbPtrBuf = dctCbPtr + i * halfWidth + j;
      twoDimDct(centerCbPtrBuf, dctCbPtrBuf);

      centerCrPtrBuf = centerCrPtr + i * halfWidth + j;
      dctCrPtrBuf = dctCrPtr + i * halfWidth + j;
      twoDimDct(centerCrPtrBuf, dctCrPtrBuf);
    }
  }

  /* Quantization */
  cv::Mat quantY(input.rows, input.cols, CV_32SC1);
  int* quantYPtr= (int*)quantY.data;
  cv::Mat quantCr(halfHeight, halfWidth, CV_32SC1);
  int* quantCrPtr= (int*)quantCr.data;
  cv::Mat quantCb(halfHeight, halfWidth, CV_32SC1);
  int* quantCbPtr= (int*)quantCb.data;
  int* quantYPtrBuf;
  int* quantCrPtrBuf;
  int* quantCbPtrBuf;

  for(int i = 0; i < height; i+=8) {
    for(int j = 0; j < width; j+=8) {
      dctYPtrBuf = dctYPtr + i * width + j;
      quantYPtrBuf = quantYPtr + i * width + j;
      applyQuantTable(dctYPtrBuf, quantYPtrBuf, quantTableY);
    }
  }
  for(int i = 0; i < halfHeight; i+=8) {
    for(int j = 0; j < halfWidth; j+=8) {
      dctCbPtrBuf = dctCbPtr + i * halfWidth + j;
      quantCbPtrBuf = quantCbPtr + i * halfWidth + j;
      applyQuantTable(dctCbPtrBuf, quantCbPtrBuf, quantTableCbCr);

      dctCrPtrBuf = dctCrPtr + i * halfWidth + j;
      quantCrPtrBuf = quantCrPtr + i * halfWidth + j;
      applyQuantTable(dctCrPtrBuf, quantCrPtrBuf, quantTableCbCr);
    }
  }

  /* Zig-Zag Scan */
  cv::Mat zigzagY(input.rows, input.cols, CV_32SC1);
  int* zigzagYPtr= (int*)zigzagY.data;
  cv::Mat zigzagCr(halfHeight, halfWidth, CV_32SC1);
  int* zigzagCrPtr= (int*)zigzagCr.data;
  cv::Mat zigzagCb(halfHeight, halfWidth, CV_32SC1);
  int* zigzagCbPtr= (int*)zigzagCb.data;
 
  /* Entropy Encoding */
  int lastYDcValue = 0;
  int lastCbDcValue = 0;
  int lastCrDcValue = 0;

#if 1
  for(int i = 0; i < halfHeight; i+=8) {
    for(int j = 0; j < halfWidth; j+=8) {
      quantYPtrBuf = quantYPtr + (i * 2) * width + (j * 2);
      code_coefficients_pass_two(quantYPtrBuf, huffCodesYDc, huffCodesSizesYDc,
          huffCodesYAc, huffCodesSizesYAc, &lastYDcValue);
      quantYPtrBuf = quantYPtr + (i * 2) * width + (j * 2) + 8;
      code_coefficients_pass_two(quantYPtrBuf, huffCodesYDc, huffCodesSizesYDc,
          huffCodesYAc, huffCodesSizesYAc, &lastYDcValue);
      quantYPtrBuf = quantYPtr + ((i * 2) + 8) * width + (j * 2);
      code_coefficients_pass_two(quantYPtrBuf, huffCodesYDc, huffCodesSizesYDc,
          huffCodesYAc, huffCodesSizesYAc, &lastYDcValue);
      quantYPtrBuf = quantYPtr + ((i * 2) + 8) * width + (j * 2) + 8;
      code_coefficients_pass_two(quantYPtrBuf, huffCodesYDc, huffCodesSizesYDc,
          huffCodesYAc, huffCodesSizesYAc, &lastYDcValue);

      quantCbPtrBuf = quantCbPtr + i * halfWidth + j;
      code_coefficients_pass_two(quantCbPtrBuf, huffCodesCbCrDc, huffCodesSizesCbCrDc,
          huffCodesCbCrAc, huffCodesSizesCbCrAc, &lastCbDcValue);
      quantCrPtrBuf = quantCrPtr + i * halfWidth + j;
      code_coefficients_pass_two(quantCrPtrBuf, huffCodesCbCrDc, huffCodesSizesCbCrDc,
          huffCodesCbCrAc, huffCodesSizesCbCrAc, &lastCrDcValue);
    }
  }
#endif

  /* Add Header */
  terminate_pass_two();

#if 1
  cv::namedWindow(windowName, CV_WINDOW_AUTOSIZE);
  cv::imshow(windowName, reverseBlockY);
  cv::waitKey(0);
#endif

#ifdef WRITE_FILES
  close();
#endif
  return 0;
}
