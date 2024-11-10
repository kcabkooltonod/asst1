作业1 - 程序1:

1.利用多个线程实现生成mandelbrot set的图片。
实现方法：
在prog1文件夹中的mandelbrotThread.cpp中实现：
```
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
```
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
```
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

