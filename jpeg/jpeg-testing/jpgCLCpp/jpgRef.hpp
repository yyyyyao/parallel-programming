#include <iostream>
#include <typeinfo>
#include <string>
#include <stdio.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "include/jpgInternal.h"

#define JPGE_OUT_BUF_SIZE 2048

class jpgRef {
  public:
    jpgRef(cv::Mat input, char* outputFileName, int quality) {
      m_input = input;
      strcpy(m_outputFileName, outputFileName);
      m_quality = quality;

      init();
      open(outputFileName);
      emit_markers();
    }

    void execute() {
      execColorConv();
      execSplitData();
      execSubSamplingH2V2();
      execBlockSplitting();
      execCenterization();
      execDct();
      execZigzag();
      execQuantization();
      execHuffman();

      terminate_pass_two();
      close();
    }

    cv::Mat getColorConvResult(void) {
      return m_colorConvNonSplit;
    }
    cv::Mat getSplitDataYResult(void) {
      return m_Y;
    }
    cv::Mat getSplitDataCbResult(void) {
      return m_Cb;
    }
    cv::Mat getSplitDataCrResult(void) {
      return m_Cr;
    }
    cv::Mat getSubSamplingCbResult(void) {
      return m_subCb;
    }
    cv::Mat getSubSamplingCrResult(void) {
      return m_subCr;;
    }
    cv::Mat getCenterizationYResult(void) {
      return m_centerY;
    }
    cv::Mat getCenterizationCbResult(void) {
      return m_centerCb;
    }
    cv::Mat getCenterizationCrResult(void) {
      return m_centerCr;
    }
    cv::Mat getDctYresult(void) {
      return m_dctY;
    }
    cv::Mat getDctCbresult(void) {
      return m_dctCb;
    }
    cv::Mat getDctCrresult(void) {
      return m_dctCr;
    }
    cv::Mat getZigzagYResult(void) {
      return m_zigzagY;
    }
    cv::Mat getZigzagCbResult(void) {
      return m_zigzagCb;
    }
    cv::Mat getZigzagCrResult(void) {
      return m_zigzagCr;
    }
    cv::Mat getQuantizationYResult(void) {
      return m_quantY;
    }
    cv::Mat getQuantizationCbResult(void) {
      return m_quantCb;
    }
    cv::Mat getQuantizationCrResult(void) {
      return m_quantCr;
    }

    void showImage(cv::Mat data) {
      cv::namedWindow("data", CV_WINDOW_AUTOSIZE);
      cv::imshow("data", data);
      cv::waitKey(0);
    }
    void execColorConv() {
      for(int i = 0; i < m_image_y; i++) {
        unsigned char* srcPtrBuf = inputPtr + i * m_image_x * 3;
        unsigned char* dstPtrBuf = colorConvNonSplitPtr + i * m_image_x * 3;
        RGB_to_YCC(dstPtrBuf, srcPtrBuf, m_image_x);
      }
    }

    void execSplitData() {
      splitToComp(colorConvNonSplitPtr,
          YPtr, CbPtr, CrPtr, m_image_x, m_image_y);
    }

    void execSubSamplingH2V2() {
      SubSamplingH2V2(CbPtr, subCbPtr, m_image_x, m_image_y);
      SubSamplingH2V2(CrPtr, subCrPtr, m_image_x, m_image_y);
    }

    void execBlockSplitting() {
      /* block Splitting Y */
      for(int i = 0; i < m_image_y; i+=8) {
        uint8* YPtrBuf = YPtr + i * m_image_x ;
        uint8* blockYPtrBuf = blockYPtr + i * m_image_x;
        blockSplitMcuLine(YPtrBuf, blockYPtrBuf, m_image_x);
      }

      /* block Splitting CbCr */
      for(int i = 0; i < m_image_half_y; i+=8) {
        uint8* subCbPtrBuf = subCbPtr + i * m_image_half_x ;
        uint8* blockCbPtrBuf = blockCbPtr + i * m_image_half_x;
        blockSplitMcuLine(subCbPtrBuf, blockCbPtrBuf, m_image_half_x);

        uint8* subCrPtrBuf = subCrPtr + i * m_image_half_x ;
        uint8* blockCrPtrBuf = blockCrPtr + i * m_image_half_x;
        blockSplitMcuLine(subCrPtrBuf, blockCrPtrBuf, m_image_half_x);
      }
    }

    void execCenterization(void) {
      centerization(blockYPtr, centerYPtr, m_image_x, m_image_y);
      centerization(blockCbPtr, centerCbPtr, m_image_half_x, m_image_half_y);
      centerization(blockCrPtr, centerCrPtr, m_image_half_x, m_image_half_y);
    }

    void execDct(void){
      for(int i = 0; i < m_image_y * m_image_x; i+=64) {
        char* centerYPtrBuf = centerYPtr + i;
        int* dctYPtrBuf = dctYPtr + i;
        twoDimDct(centerYPtrBuf, dctYPtrBuf);

#ifdef PRINT_BLOCK
        printf("------centerization result Y%d-------\n", i / 64);
        showBlocks(centerYPtrBuf);
        printf("------dct result Y%d-------\n", i / 64);
        showBlocks(dctYPtrBuf);
#endif
      }

      for(int i = 0; i < m_image_half_y * m_image_half_x; i+=64) {
        char* centerCbPtrBuf = centerCbPtr + i;
        int* dctCbPtrBuf = dctCbPtr + i;
        twoDimDct(centerCbPtrBuf, dctCbPtrBuf);

        char* centerCrPtrBuf = centerCrPtr + i;
        int* dctCrPtrBuf = dctCrPtr + i;
        twoDimDct(centerCrPtrBuf, dctCrPtrBuf);
      }
    }

    void execZigzag(void) {
      for(int i = 0; i < m_image_y * m_image_x; i+=64) {
        int* dctYPtrBuf = dctYPtr + i;
        int* zigzagYPtrBuf = zigzagYPtr + i;
        zigZagScan(dctYPtrBuf, zigzagYPtrBuf);
      }

      for(int i = 0; i < m_image_half_y * m_image_half_x; i+=64) {
        int* dctCbPtrBuf = dctCbPtr + i;
        int* zigzagCbPtrBuf = zigzagCbPtr + i;
        zigZagScan(dctCbPtrBuf, zigzagCbPtrBuf);
        int* dctCrPtrBuf = dctCrPtr + i;
        int* zigzagCrPtrBuf = zigzagCrPtr + i;
        zigZagScan(dctCrPtrBuf, zigzagCrPtrBuf);
      }
    }

    void execQuantization(void) {
      for(int i = 0; i < m_image_y * m_image_x; i+=64) {
        int* zigzagYPtrBuf = zigzagYPtr + i;
        int* quantYPtrBuf = quantYPtr + i;
        applyQuantTable(zigzagYPtrBuf, quantYPtrBuf, m_quantTableY);

#ifdef PRINT_BLOCK
        printf("-----quant & zizag Y:%d--------\n", i / 64);
        showBlocks(quantYPtrBuf);  
#endif
      }

      for(int i = 0; i < m_image_half_y * m_image_half_x; i+=64) {
        int* zigzagCbPtrBuf = zigzagCbPtr + i;
        int* quantCbPtrBuf = quantCbPtr + i;
        applyQuantTable(zigzagCbPtrBuf, quantCbPtrBuf, m_quantTableCbCr);

        int* zigzagCrPtrBuf = zigzagCrPtr + i;
        int* quantCrPtrBuf = quantCrPtr + i;
        applyQuantTable(zigzagCrPtrBuf, quantCrPtrBuf, m_quantTableCbCr);

#ifdef PRINT_BLOCK
        printf("-----quant & zizag Cb:%d--------\n", i / 64);
        showBlocks(quantCbPtrBuf);
        printf("-----quant & zizag Cr:%d--------\n", i / 64);
        showBlocks(quantCbPtrBuf);
#endif
      }

    }

    void execHuffman(void) {
      for(int i = 0; i < m_mcuY; i++) {
        for(int j = 0; j < m_mcuX; j++) {
          int* quantYPtrBuf = quantYPtr + (i *2) * m_lineOfData * 2 + j * 128;
          code_coefficients_pass_two(quantYPtrBuf, m_huffCodesYDc, m_huffCodesSizesYDc,
              m_huffCodesYAc, m_huffCodesSizesYAc, &m_lastYDcValue);
          quantYPtrBuf = quantYPtr + (i * 2) * m_lineOfData * 2 + j * 128 + 64;
          code_coefficients_pass_two(quantYPtrBuf, m_huffCodesYDc, m_huffCodesSizesYDc,
              m_huffCodesYAc, m_huffCodesSizesYAc, &m_lastYDcValue);
          quantYPtrBuf = quantYPtr + (i * 2 + 1) * m_lineOfData * 2 + j * 128;
          code_coefficients_pass_two(quantYPtrBuf, m_huffCodesYDc, m_huffCodesSizesYDc,
              m_huffCodesYAc, m_huffCodesSizesYAc, &m_lastYDcValue);

          quantYPtrBuf = quantYPtr + (i * 2 + 1) * m_lineOfData * 2 + j * 128 + 64;
          code_coefficients_pass_two(quantYPtrBuf, m_huffCodesYDc, m_huffCodesSizesYDc,
              m_huffCodesYAc, m_huffCodesSizesYAc, &m_lastYDcValue);

          int* quantCbPtrBuf = quantCbPtr + i * m_lineOfData + j * 64;
          code_coefficients_pass_two(quantCbPtrBuf, m_huffCodesCbCrDc, m_huffCodesSizesCbCrDc,
              m_huffCodesCbCrAc, m_huffCodesSizesCbCrAc, &m_lastCbDcValue);

          int* quantCrPtrBuf = quantCrPtr + i * m_lineOfData + j * 64;
          code_coefficients_pass_two(quantCrPtrBuf, m_huffCodesCbCrDc, m_huffCodesSizesCbCrDc,
              m_huffCodesCbCrAc, m_huffCodesSizesCbCrAc, &m_lastCrDcValue);
        }
      }
    }
    
  private:
    /* general datas */
    char m_outputFileName[100];
    int m_image_x, m_image_y;
    int m_image_half_x, m_image_half_y;
    int m_quality;
    int m_mcuX;
    int m_mcuY;
    int m_lineOfData;

    cv::Mat m_input;
    uint8 *inputPtr;
    cv::Mat m_colorConvNonSplit;
    uint8 *colorConvNonSplitPtr;
    cv::Mat m_Y, m_Cb, m_Cr;
    uint8 *YPtr, *CbPtr, *CrPtr;
    cv::Mat m_subCb, m_subCr;
    uint8 *subCbPtr, *subCrPtr;
    cv::Mat m_blockY, m_blockCb, m_blockCr;
    uint8 *blockYPtr, *blockCbPtr, *blockCrPtr;
    cv::Mat m_centerY, m_centerCb, m_centerCr;
    int8 *centerYPtr, *centerCbPtr, *centerCrPtr;
    cv::Mat m_dctY, m_dctCb, m_dctCr;
    int *dctYPtr, *dctCbPtr, *dctCrPtr;
    cv::Mat m_zigzagY,  m_zigzagCb, m_zigzagCr;
    int *zigzagYPtr, *zigzagCbPtr, *zigzagCrPtr;
    cv::Mat m_quantY, m_quantCb, m_quantCr;
    int *quantYPtr, *quantCbPtr, *quantCrPtr;

    int* m_quantization_tables[3];
    int m_quantTableY[64];
    int m_quantTableCbCr[64];

    unsigned int m_huffCodesYDc[256];
    unsigned int m_huffCodesYAc[256];
    unsigned int m_huffCodesCbCrDc[256];
    unsigned int m_huffCodesCbCrAc[256];
    uint8 m_huffCodesSizesYDc[256];
    uint8 m_huffCodesSizesYAc[256];
    uint8 m_huffCodesSizesCbCrDc[256];
    uint8 m_huffCodesSizesCbCrAc[256];

    int m_lastYDcValue;
    int m_lastCbDcValue;
    int m_lastCrDcValue;

    /* for writing bits variables */
    unsigned int m_out_buf_left;
    int m_bits_in;
    uint8 m_out_buf[JPGE_OUT_BUF_SIZE];
    //uint8* m_pOut_buf = m_out_buf;
    uint8* m_pOut_buf;
    bool m_all_stream_writes_succeeded;
    unsigned int m_bit_buffer;
    bool m_bStatus;
    FILE* m_pFile;

    void init() { /* just init the member variables. */
      m_out_buf_left = JPGE_OUT_BUF_SIZE;
      memset(m_out_buf, 0, JPGE_OUT_BUF_SIZE);
      m_pFile = NULL;
      m_pOut_buf = m_out_buf;
      m_bits_in = 0;
      m_bit_buffer = 0;
      m_bStatus = true;
      m_image_x = m_input.cols;
      m_image_y = m_input.rows;
      m_image_half_x = m_image_x / 2;
      m_image_half_y = m_image_y / 2;
      m_quantization_tables[0] = m_quantTableY;
      m_quantization_tables[1] = m_quantTableCbCr;
      m_quantization_tables[2] = m_quantTableCbCr;

      m_lastYDcValue = 0;
      m_lastCbDcValue = 0;
      m_lastCrDcValue = 0;
      m_mcuX = m_image_x / 16;
      m_mcuY = m_image_y / 16;
      m_lineOfData = 64 * m_mcuX;

      matInit();

      genTables();
    }

    void genTables() {
      /* compute quantization tables. */
      calcQuantTableByQuality(s_std_lum_quant, m_quantTableY, m_quality);
      calcQuantTableByQuality(s_std_croma_quant, m_quantTableCbCr, m_quality);
      /* compute huffman Y Table. */
      compute_huffman_table(m_huffCodesYDc, m_huffCodesSizesYDc, s_dc_lum_bits, s_dc_lum_val);
      compute_huffman_table(m_huffCodesYAc, m_huffCodesSizesYAc, s_ac_lum_bits, s_ac_lum_val);
      /* compute huffman CbCr Table. */
      compute_huffman_table(m_huffCodesCbCrDc, m_huffCodesSizesCbCrDc, s_dc_chroma_bits, s_dc_chroma_val);
      compute_huffman_table(m_huffCodesCbCrAc, m_huffCodesSizesCbCrAc, s_ac_chroma_bits, s_ac_chroma_val);

    }

    void matInit() {
      inputPtr = (uint8*)m_input.data;
      m_colorConvNonSplit.create(m_image_x, m_image_y, CV_8UC3);
      colorConvNonSplitPtr = (uint8*)m_colorConvNonSplit.data;
      m_Y.create(m_image_x, m_image_y, CV_8UC1);
      YPtr = (uint8*)m_Y.data;
      m_Cb.create(m_image_x, m_image_y, CV_8UC1);
      CbPtr = (uint8*)m_Cb.data;
      m_Cr.create(m_image_x, m_image_y, CV_8UC1);
      CrPtr = (uint8*)m_Cr.data;

      m_subCb.create(m_image_half_x, m_image_half_y, CV_8UC1);
      subCbPtr = (uint8*)m_subCb.data;
      m_subCr.create(m_image_half_x, m_image_half_y, CV_8UC1);
      subCrPtr = (uint8*)m_subCr.data;

      m_blockY.create(m_image_x, m_image_y, CV_8UC1);
      blockYPtr = (uint8*)m_blockY.data;
      m_blockCb.create(m_image_half_x, m_image_half_y, CV_8UC1);
      blockCbPtr = (uint8*)m_blockCb.data;
      m_blockCr.create(m_image_half_x, m_image_half_y, CV_8UC1);
      blockCrPtr = (uint8*)m_blockCr.data;

      m_centerY.create(m_image_x, m_image_y, CV_8UC1);
      centerYPtr = (int8*)m_centerY.data;
      m_centerCb.create(m_image_half_x, m_image_half_y, CV_8UC1);
      centerCbPtr = (int8*)m_centerCb.data;
      m_centerCr.create(m_image_half_x, m_image_half_y, CV_8UC1);
      centerCrPtr = (int8*)m_centerCr.data;

      m_dctY.create(m_image_x, m_image_y, CV_32SC1);
      dctYPtr = (int*)m_dctY.data;
      m_dctCb.create(m_image_half_x, m_image_half_y, CV_32SC1);
      dctCbPtr = (int*)m_dctCb.data;
      m_dctCr.create(m_image_half_x, m_image_half_y, CV_32SC1);
      dctCrPtr = (int*)m_dctCr.data;

      m_zigzagY.create(m_image_x, m_image_y, CV_32SC1);
      zigzagYPtr = (int*)m_zigzagY.data;
      m_zigzagCb.create(m_image_half_x, m_image_half_y, CV_32SC1);
      zigzagCbPtr = (int*)m_zigzagCb.data;
      m_zigzagCr.create(m_image_half_x, m_image_half_y, CV_32SC1);
      zigzagCrPtr = (int*)m_zigzagCr.data;

      m_quantY.create(m_image_x, m_image_y, CV_32SC1);
      quantYPtr = (int*)m_quantY.data;
      m_quantCb.create(m_image_half_x, m_image_half_y, CV_32SC1);
      quantCbPtr = (int*)m_quantCb.data;
      m_quantCr.create(m_image_half_x, m_image_half_y, CV_32SC1);
      quantCrPtr = (int*)m_quantCr.data;
    }

    template<class T> void showBlocks(T b) {
      for(int i = 0; i < 8; i++) {
        for(int j = 0; j < 8; j++) {
          printf("%4d", b[i * 8 + j]);
        }
        printf("\n");
      }
      return;
    }

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

    bool put_buf(void* pBuf, int len)
    {
      bool temp = false;
      temp = (fwrite(pBuf, len, 1, m_pFile) == 1);
      if(temp == false) {
        printf("error!@put_buf\n");
        exit(1);
      }
      m_bStatus = m_bStatus && temp;
      return m_bStatus;
    }

    template<class T> inline bool put_obj(T& obj) { return put_buf(&obj, sizeof(T)); }


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

#if 0
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
      while (m_bits_in >= 8)/* if buffer length exceeds 8-bit. */
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
      put_bits(0x7F, 7); /* if some bits remain, weshould put away. blank bits should be 1. sor we use 0x7F, 7*/
      flush_output_buffer(); /* after the put_bits we have to flush. */
      emit_marker(M_EOI); /* emit_marker do flush in this function. */
      return true;
    }

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
      }
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

#if 0
    void applyQuantTable(const int* pSrc, int* pDst, int* quantTable)
    {
      for(int i = 0; i < 64; i++) {
        pDst[i] = pSrc[i] / quantTable[i];
      }
    }
#else
    void applyQuantTable(const int* pSrc, int* pDst, int* quantTable)
    {
      int *q = quantTable; 
      for (int i = 0; i < 64; i++)
      {
        int j = pSrc[i]; /* apply quantization. and zig-zag. */
        //int j = pSrc[s_zag[i]]; /* apply quantization. and zig-zag. */
        if (j < 0) /* if DCT data is minus. */
        {
          if ((j = -j + (*q >> 1)) < *q)
            *pDst++ = 0;
          else
            *pDst++ = static_cast<int>(-(j / *q));
        }
        else
        {
          if ((j = j + (*q >> 1)) < *q)
            *pDst++ = 0;
          else
            *pDst++ = static_cast<int>((j / *q));
        }
        q++;
      }
    }
#endif

#if 0
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
#else

#define JPGE_MAX(a,b) (((a)>(b))?(a):(b))                  
#define JPGE_MIN(a,b) (((a)<(b))?(a):(b))
    // Quantization table generation.
    void calcQuantTableByQuality(int *pSrc, int* pDst, int quality)
    {
      int32 q;
      if (quality < 50)
        q = 5000 / quality;
      else
        q = 200 - quality * 2;
      for (int i = 0; i < 64; i++)
      {
        int32 j = *pSrc++; j = (j * q + 50L) / 100L;
        *pDst++ = JPGE_MIN(JPGE_MAX(j, 1), 255);
      }
    }
#endif

    void centerization(const uint8* pSrc, int8* pDst, int width, int height)
    {
      for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
          pDst[i * width + j] = (int8)(pSrc[i * width + j] - 128);
        }
      }

    }

#if 0
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
#else
    // Forward DCT - DCT derived from jfdctint.
    enum { CONST_BITS = 13, ROW_BITS = 2 };
#define DCT_DESCALE(x, n) (((x) + (((int32)1) << ((n) - 1))) >> (n))
#define DCT_MUL(var, c) ((var) * (c))
#define DCT1D(s0, s1, s2, s3, s4, s5, s6, s7) \
      int32 t0 = s0 + s7, t7 = s0 - s7, t1 = s1 + s6, t6 = s1 - s6, t2 = s2 + s5, t5 = s2 - s5, t3 = s3 + s4, t4 = s3 - s4; \
      int32 t10 = t0 + t3, t13 = t0 - t3, t11 = t1 + t2, t12 = t1 - t2; \
      int32 u1 = DCT_MUL(t12 + t13, 4433); \
      s2 = u1 + DCT_MUL(t13, 6270); \
      s6 = u1 + DCT_MUL(t12, -15137); \
      u1 = t4 + t7; \
      int32 u2 = t5 + t6, u3 = t4 + t6, u4 = t5 + t7; \
      int32 z5 = DCT_MUL(u3 + u4, 9633); \
      t4 = DCT_MUL(t4, 2446); t5 = DCT_MUL(t5, 16819); \
      t6 = DCT_MUL(t6, 25172); t7 = DCT_MUL(t7, 12299); \
      u1 = DCT_MUL(u1, -7373); u2 = DCT_MUL(u2, -20995); \
      u3 = DCT_MUL(u3, -16069); u4 = DCT_MUL(u4, -3196); \
      u3 += z5; u4 += z5; \
      s0 = t10 + t11; s1 = t7 + u1 + u4; s3 = t6 + u2 + u3; s4 = t10 - t11; s5 = t5 + u2 + u4; s7 = t4 + u1 + u3;

    void castInt8_to_Int32(int8* pSrc, int* pDst) {
      for(int i = 0; i < 64; i++) {
        pDst[i] = (int)pSrc[i];
      }
    }
    void twoDimDct(int8* pSrc, int* pDst)
    {
      int castSrc[64];
      castInt8_to_Int32(pSrc, castSrc);

      int* p = castSrc;
      int c, *q = p;
      for (c = 7; c >= 0; c--, q += 8)
      {
        int s0 = q[0], s1 = q[1], s2 = q[2], s3 = q[3], s4 = q[4], s5 = q[5], s6 = q[6], s7 = q[7];
        DCT1D(s0, s1, s2, s3, s4, s5, s6, s7);
        q[0] = s0 << ROW_BITS; q[1] = DCT_DESCALE(s1, CONST_BITS-ROW_BITS); q[2] = DCT_DESCALE(s2, CONST_BITS-ROW_BITS); q[3] = DCT_DESCALE(s3, CONST_BITS-ROW_BITS);
        q[4] = s4 << ROW_BITS; q[5] = DCT_DESCALE(s5, CONST_BITS-ROW_BITS); q[6] = DCT_DESCALE(s6, CONST_BITS-ROW_BITS); q[7] = DCT_DESCALE(s7, CONST_BITS-ROW_BITS);
      }
      for (q = p, c = 7; c >= 0; c--, q++)
      {
        int s0 = q[0*8], s1 = q[1*8], s2 = q[2*8], s3 = q[3*8], s4 = q[4*8], s5 = q[5*8], s6 = q[6*8], s7 = q[7*8];
        DCT1D(s0, s1, s2, s3, s4, s5, s6, s7);
        q[0*8] = DCT_DESCALE(s0, ROW_BITS+3); q[1*8] = DCT_DESCALE(s1, CONST_BITS+ROW_BITS+3); q[2*8] = DCT_DESCALE(s2, CONST_BITS+ROW_BITS+3); q[3*8] = DCT_DESCALE(s3, CONST_BITS+ROW_BITS+3);
        q[4*8] = DCT_DESCALE(s4, ROW_BITS+3); q[5*8] = DCT_DESCALE(s5, CONST_BITS+ROW_BITS+3); q[6*8] = DCT_DESCALE(s6, CONST_BITS+ROW_BITS+3); q[7*8] = DCT_DESCALE(s7, CONST_BITS+ROW_BITS+3);
      }

      for(int i = 0; i < 64; i++) {
        pDst[i] = castSrc[i];
      }
    }
#endif

    void reverseBlockSplit(const uint8* pSrc, uint8* pDst, int width, int height)
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

    static inline uint8 clamp(int i) { if (static_cast<uint>(i) > 255U) { if (i < 0) i = 0; else if (i > 255) i = 255; } return static_cast<uint8>(i); }
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
};
