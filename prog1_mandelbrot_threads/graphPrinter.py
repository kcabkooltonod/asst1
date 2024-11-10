import matplotlib.pyplot as plt

# 读取数据
threads = []
speedup = []

with open('speedup_data.txt', 'r') as file:
    next(file)  # Skip header
    for line in file:
        t, s = map(float, line.split())
        threads.append(int(t))
        speedup.append(s)

# 绘制折线图
plt.figure(figsize=(10, 6))
plt.plot(threads, speedup, marker='o', color='b', linestyle='-')
plt.xlabel("Number of Threads")
plt.ylabel("Speedup (x)")
plt.title("Speedup vs Number of Threads")
plt.grid(True)
plt.xticks(threads)

# 保存图片为 speedup_plot.png
plt.savefig("speedup_plot1.png")