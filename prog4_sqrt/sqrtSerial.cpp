#include <immintrin.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

void sqrtSerial(int N,
                float initialGuess,
                float values[],
                float output[])
{

    static const float kThreshold = 0.00001f;

    for (int i=0; i<N; i++) {

        float x = values[i];
        float guess = initialGuess;

        float error = fabs(guess * guess * x - 1.f);

        while (error > kThreshold) {
            guess = (3.f * guess - x * guess * guess * guess) * 0.5f;
            error = fabs(guess * guess * x - 1.f);
        }

        output[i] = x * guess;
    }
}

void sqrtAVX2(int N,
              float initialGuess,
              float values[],
              float output[])
{
    static const float kThreshold = 0.00001f;
    __m256 threshold = _mm256_set1_ps(kThreshold);
    __m256 three = _mm256_set1_ps(3.0f);
    __m256 half = _mm256_set1_ps(0.5f);

    int i = 0;
    for (; i <= N - 8; i += 8) {
        __m256 x = _mm256_loadu_ps(&values[i]);
        __m256 guess = _mm256_set1_ps(initialGuess);
        
        __m256 error = _mm256_sub_ps(_mm256_mul_ps(_mm256_mul_ps(guess, guess), x), _mm256_set1_ps(1.0f));
        // 首先，生成一个带有符号的-0.f，符号位是1其余为0
        // 然后就是牛逼的andnot指令，先对前一个寄存器取反，然后与第二个寄存器做与运算
        // 由于取反后的-0.f的符号位为0，其余为一，所以保留了error的数值并且删除了符号位，达到绝对值的效果
        error = _mm256_andnot_ps(_mm256_set1_ps(-0.0f), error);
        
        //movemask指令提取8个数字的符号位并返回一个8bit，只有当全0的时候才跳出循环，也就是所以
        while (_mm256_movemask_ps(_mm256_cmp_ps(error, threshold, _CMP_GT_OQ))) {
            guess = _mm256_mul_ps(_mm256_sub_ps(_mm256_mul_ps(three, guess), 
                            _mm256_mul_ps(_mm256_mul_ps(x, guess), _mm256_mul_ps(guess, guess))), half);
            
            error = _mm256_sub_ps(_mm256_mul_ps(_mm256_mul_ps(guess, guess), x), _mm256_set1_ps(1.0f));
            error = _mm256_andnot_ps(_mm256_set1_ps(-0.0f), error);  // fabs
        }
        
        __m256 result = _mm256_mul_ps(x, guess);
        _mm256_storeu_ps(&output[i], result);
    }

    // Handle remaining elements
    for (; i < N; i++) {
        float x = values[i];
        float guess = initialGuess;

        float error = fabs(guess * guess * x - 1.f);
        while (error > kThreshold) {
            guess = (3.f * guess - x * guess * guess * guess) * 0.5f;
            error = fabs(guess * guess * x - 1.f);
        }
        output[i] = x * guess;
    }
}