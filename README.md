## 作业1 - 程序1:

1.利用多个线程实现生成mandelbrot set的图片。
实现方法：
在prog1文件夹中的mandelbrotThread.cpp中实现：
```C++
int rowsPerThread = args->height / args->numThreads;
    int startRow = args->threadId * rowsPerThread;
    int numRows = rowsPerThread;
    if (args->threadId + 1 == args->numThreads) {
        numRows = args->height - startRow;
    }
    mandelbrotSerial(args->x0, args->y0, args->x1, args->y1, args->width, args->height, 
                        startRow, numRows, args->maxIterations, args->output);    
```

2.更改线程个数，并统计时间后会发现加速比并不等于线程的个数，甚至使用3个线程比使用2个线程还要慢，这是因为上面的方法是从上到下依次均匀分块，由于需要生成的点在图像上并不均匀而导致了线程分配到的工作是不均衡的，从而分到中间部分的线程会工作更久。


3.在调用
`mandelbrotSerial()`前后加上统计时间的函数后会发现中间的线程明显工作时间更久。
```C++
int rowsPerThread = args->height / args->numThreads;
    int startRow = args->threadId * rowsPerThread;
    int numRows = rowsPerThread;
    if (args->threadId + 1 == args->numThreads) {
        numRows = args->height - startRow;
    }
    // 加入计时
    double startTime = CycleTimer::currentSeconds();
    mandelbrotSerial(args->x0, args->y0, args->x1, args->y1, args->width, args->height, 
                        startRow, numRows, args->maxIterations, args->output);    
    double endTime = CycleTimer::currentSeconds();

    printf("[mandelbrot thread %d]:\t\t[%.3f] ms\n", args->threadId, (endTime - startTime) * 1000);
```

4.解决办法：放弃之前的连续分配方法，让每一个线程分到的行数根据线程的个数进行跨步。
```C++
// 新加的分配方法
void mandelbrotSerialNew(
    float x0, float y0, float x1, float y1,
    int width, int height,
    int startRow, int step,
    int maxIterations,
    int output[])
{
    float dx = (x1 - x0) / width;
    float dy = (y1 - y0) / height;

    // int endRow = startRow + totalRows;

    for (int j = startRow; j < height; j += step) {
        for (int i = 0; i < width; ++i) {
            float x = x0 + i * dx;
            float y = y0 + j * dy;

            int index = (j * width + i);
            output[index] = mandel(x, y, maxIterations);
        }
    }
}

// 更改后的workerThreadStart()方法
    double startTime = CycleTimer::currentSeconds();
    mandelbrotSerialNew(args->x0, args->y0, args->x1, args->y1, args->width, args->height, 
                        args->threadId, args->numThreads, args->maxIterations, args->output);    
    double endTime = CycleTimer::currentSeconds();
    printf("[mandelbrot thread %d]:\t\t[%.3f] ms\n", args->threadId, (endTime - startTime) * 1000);

```

5.把线程加到16后，我的机器仍然可以继续增加，最多增加到了12倍左右，这是由机器的核心数量和每个核心可以交错执行线程的数量决定的。比如我的机器是10个core、每个core上的hyper-thread是2.
可以用`lscpu`命令去查看cpu的信息。

## 程序二

1.用向量指令SIMD实现```clampedExpVector()```方法，并实现边角料的处理。
```C++
void clampedExpVector(float* values, int* exponents, float* output, int N) {

  //
  // CS149 STUDENTS TODO: Implement your vectorized version of
  // clampedExpSerial() here.
  //
  // Your solution should work for any value of
  // N and VECTOR_WIDTH, not just when VECTOR_WIDTH divides N
  //
  int elemLeft = N % VECTOR_WIDTH;
  __cs149_vec_float value;
  __cs149_vec_int exponent;
  __cs149_vec_float result;
  __cs149_vec_int zero = _cs149_vset_int(0);
  __cs149_vec_int one = _cs149_vset_int(1);
  __cs149_vec_float max = _cs149_vset_float(9.999999f);

  __cs149_mask maskAll, maskIsNotZero, maskIsMax;

  maskAll = _cs149_init_ones();

  for(int i = 0; i < N - elemLeft; i += VECTOR_WIDTH) {
    result = _cs149_vset_float(1.f);

    _cs149_vload_float(value, values+i, maskAll);
    _cs149_vload_int(exponent, exponents+i, maskAll);

    _cs149_vgt_int(maskIsNotZero, exponent, zero, maskAll);

    while (_cs149_cntbits(maskIsNotZero)) {
      _cs149_vmult_float(result, result, value, maskIsNotZero);
      _cs149_vsub_int(exponent, exponent, one, maskIsNotZero);
      _cs149_vgt_int(maskIsNotZero, exponent, zero, maskIsNotZero);
    }

    _cs149_vgt_float(maskIsMax, result, max, maskAll);
    _cs149_vset_float(result, 9.999999f, maskIsMax);

    _cs149_vstore_float(output+i, result, maskAll);

  }

  if (elemLeft) {
    clampedExpSerial(values + (N - elemLeft), exponents + (N - elemLeft), output + (N - elemLeft), elemLeft);
  }
  
}
```

2.通过改变SIMD指令中的向量长度，会发现随着向量长度的增加，vector utilization使用效率会减少，虽然不是很多，但是还是在稳定减少的，说明了随着向量长度的增加，一个向量中就会存在更多的分化现象，从而导致not coherent control flow的概率越来越大。


3.利用向量指令实现数组求和，提示使用hadd方法和interleave方法。hadd方法的实现：C++
```
void _cs149_hadd(__cs149_vec<T> &vecResult, __cs149_vec<T> &vec) {
  for (int i=0; i<VECTOR_WIDTH/2; i++) {
    T result = vec.value[2*i] + vec.value[2*i+1];
    vecResult.value[2 * i] = result;
    vecResult.value[2 * i + 1] = result;
  }
}
```

interleave方法的实现：
```C++
void _cs149_interleave(__cs149_vec<T> &vecResult, __cs149_vec<T> &vec) {
  for (int i=0; i<VECTOR_WIDTH; i++) {
    int index = i < VECTOR_WIDTH/2 ? (2 * i) : (2 * (i - VECTOR_WIDTH/2) + 1);
    vecResult.value[i] = vec.value[index];
  }
}
```
可以大致看出hadd方法是将向量中相邻的两个数据替换为它们的和，例如一个长度为4的向量值分别为：1，2，3，4经过hadd后会变成3，3，7，7

interleave方法就是将hadd中不同的数据全放到前半部分去，以便下一次进行循环。相当于每一次经过这两个操作后，需要求和的数组长度就会变成前一次的一半。

因此，就有了以下```arraySumVector()```方法的实现
```C++
float arraySumVector(float* values, int N) {
  
  //
  // CS149 STUDENTS TODO: Implement your vectorized version of arraySumSerial here
  //
  __cs149_vec_float sum = _cs149_vset_float(0.f);
  __cs149_vec_float temp;
  __cs149_mask maskAll;

  maskAll = _cs149_init_ones();

  for (int i = 0; i < N; i += VECTOR_WIDTH) {
    _cs149_vload_float(temp, values + i, maskAll);
    _cs149_vadd_float(sum, sum, temp, maskAll);
  }

  for (int i = VECTOR_WIDTH; i > 1; i /= 2) {
    _cs149_hadd_float(sum, sum);
    _cs149_interleave_float(sum, sum);
  }

  float output[VECTOR_WIDTH];
  _cs149_vstore_float(output, sum, maskAll);

  return output[0];

}
```

## 程序三

1.1.在安装好ispc后可以进入到prog3的文件夹里进行make了。然后运行可执行文件后会发现提升的速度只有一半左右也就是4x，但明明SIMD中的向量长度可是8个浮点数。这跟prog2中的情况是一样的，由于生成的图片中黑色和白色的点并没有完全连续的分布在一起，而SIMD指令最好是能处理连续的相同操作，即避免分支的情况。因此与理想的加速比相差甚远。

2.1.使用--task选项后会发现加速比提升了一点（2个task的情况），task其实是在thread之上的一个对象，是一种轻量级的逻辑并行对象，一个thread上可以被分配多个tasks并行执行。随着task的数量增加，加速比会飞快的提升，我的机器上最多提升了32x。

2.2.分配的task增加后，当到达40的时候取得了最大加速32x。不同机器可能不一样。

2.3.task可以是无限多的，但也不能太多，因为太多task把要生成的图片区域划分的是在太小的时候会出现错误：
```
Wrote image file mandelbrot-task-ispc.ppm
Mismatch : [0][108], Expected : 1, Actual : 0
Error : ISPC output differs from sequential output
```
可能是分太细了吧。

thread的数量是有限制的，我的机器上是32个。

