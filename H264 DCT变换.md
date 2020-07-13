# `H264` `DCT`变换

## 基本概念

变换的是用来减少图像编码的动态范围，将图像的时域信号变换为频域信号。频域中，图像的能量大部分集中在低频区域，相对于时域信号，码率有所下降。

`H264`中对图像或残差采用的是$4\times4$**整数**离散余弦变换技术。`H264`编码将变换编码和量化两个过程中的乘法合二为一，并进一步的采用整数的运算形式，减少编解码过程中的计算量，提高图像压缩的实时性。整数变化和量化的过程如下所示：

![image-20200620145516980](/Users/mac/Library/Application Support/typora-user-images/image-20200620145516980.png)

从流程图中可以看出第一次量化缩放之后会判断输入的$4 \times 4$块是否为色度块或者采用的是帧内$16 \times16$的预测模式，如果是的话就会在进行$HDM$变换，然后再次量化（对于这一块放在量化的模块在讨论）。

## 一维$DCT$变换

基本公式：
$$
y_k=C_k\sum_{n=0}^{N-1}x_n\cos\frac{(2n+1)k\pi}{2N}
$$

$$
C_k=\begin{cases} \sqrt\frac{1}{N}, k=0\\ \sqrt\frac{2}{N}, k=1,2,...N-1\end{cases}
$$

每个$DCT$系数$y_k$确定信号$x_n$在相应频率点上的贡献。最低的系数（$k=0$）为$DC$系数，代表信号的平均值，也称为直流分量，其它的系数称为交流系数$AC$。

一维$N$点离散余弦逆变换($IDCT$)：
$$
x_n=\sum_{k=0}^{N-1}C_ky_k\cos\frac{(2n+1)k\pi}{2N}
$$
来直接上一维$DCT$变换的代码：

```c
#define PI 3.1415926
#define C0 ((double)1/sqrt((double)len))
#define C1 ((sqrt((double)2))/sqrt((double)len))

void oneDct(double* pDct, uint8_t* pData, int32_t len)
{
    int i;
    int j;
    int val = 0;
    for (j = 0; j < len; j++) {
        val += pData[j] * cos(0);
    }
    pDct[0] = (double)C0 * (double)val;
    for (i = 1; i < len; i++) {
        val = 0;
        for (j = 0; j < len; j++) {
            val += pData[j] * cos((2 * j + 1) * i * PI / (2 * len));
        }
        pDct[i] = (double)C1 * (double)val;
    }
}
```

一维$IDCT$变换代码：

```c
void oneIDct(double* pDct, double* pData, int32_t len)
{
    int i, j;
    for (i = 0; i < len; i++) {
        double val = 0;
        val += C0 * pDct[0] * cos(0);
        for (int j = 1; j < len; j++) {
            val += C1 * pDct[j] * cos((2 * i + 1) * j * PI / (2 * len));
        }
        pData[i] = val;
    }
}
```



## 二维$DCT$变换

二维$N\times N$图像的$DCT$变换可以理解为先对图像块的每行进行一维$DCT$，然后对经过行变换的块的列在应用一维$DCT$变换。数学公式：
$$
Y_{mn}=C_mC_n\sum_{i=0}^{N-1}\sum_{j=0}^{N-1}X_{ij}\cos\frac{(2j+1)n\pi}{2N}\cos\frac{(2i+1)m\pi}{2N}
$$
其中$X_{ij}$是图像块$X$中$i$行$j$列图像或者是残差值，$Y_{mn}$是变换结果矩阵中$Y$相应频率点上的$DCT$系数。可以使用如下的矩阵形式表示：
$$
Y=AXA^T
$$
二维$IDCT$公式：
$$
X_{ij}=\sum_{i=0}^{N-1}\sum_{j=0}^{N-1}C_mC_n\cos\frac{(2j+1)n\pi}{2N}\cos\frac{(2i+1)m\pi}{2N}
$$

$$
X=A^TYA
$$

其中$N\times N$的变换矩阵$A$中的系数：
$$
A_{ij}=C_i\cos\frac{(2j+1)i\pi}{2N}
$$
二维$DCT$变换的$C$代码：

```c++
void DCT(const unsigned char block[N][N], double DCT_block[N][N]){

    // 申请临时空间
    double *tmp = new double[N*N];
    double *coff = new double[N*N];
    memset(tmp, 0, sizeof(double)*N*N);
    // 设置DCT系数
    coff[0] = 1.0 / sqrt((double)N);
    for (int m = 1; m < N; m++){
        coff[m] = sqrt((double)2) / sqrt((double)N);
    }

    // 实现DCT变换
    for (int m = 0; m < N; m++){
        for (int l = 0; l < N; l++){
            for (int x = 0; x < N; x++){
                tmp[m*N + l] += coff[l] * block[m][x] * cos((2 * x + 1)*PI*l / (2 * N));
            }
        }
    }
    for (int k = 0; k < N; k++){
        for (int l = 0; l < N; l++){
            for (int x = 0; x < N; x++){
                DCT_block[k][l] += coff[k] * tmp[x*N + l] * cos((2 * x + 1)*PI*k / (2 * N));
            }
        }
    }

    // 销毁分配的空间
    delete[]tmp;
    delete[]coff;
}
```

二维$IDCT$变换的$C$代码：

```c++
void IDCT(const double block[N][N], double IDCT_block [N][N]){

    // 申请临时空间
    double *tmp = new double[N*N];
    double *coff = new double[N*N];
    memset(tmp, 0, sizeof(double)*N*N);
    // 设置IDCT系数
    coff[0] = 1.0 / sqrt((double)N);
    for (int m = 1; m < N; m++){
        coff[m] = sqrt((double)2) / sqrt((double)N);
    }

    // 实现IDCT变换
    for (int k = 0; k < N; k++){
        for (int n = 0; n < N; n++){
            for (int x = 0; x < N; x++){
                tmp[k*N + n] += coff[x] * block[k][x] * cos((2 * n + 1)*x*PI / 2 / N);
            }
        }
    }
    for (int m = 0; m < N; m++){
        for (int n = 0; n < N; n++){
            for (int x = 0; x < N; x++){
                IDCT_block[m][n] += coff[x] * tmp[x*N + n] * cos((2 * m + 1)*x*PI / 2 / N);
            }
        }
    }

    // 销毁分配的空间
    delete[]tmp;
    delete[]coff;
}
```

## $H264$中二维$DCT$变换

对于$H264$中$4\times4$的图像块进行操作，对应的$4\times4$变换矩阵$A$的系数：

<img src="/Users/mac/Library/Application Support/typora-user-images/image-20200621142305982.png" alt="image-20200621142305982" style="zoom: 33%;" />

设$a=\frac{1}{2}$，$b=\sqrt{\frac{1}{2}}\cos(\frac{\pi}{8})$及$c=\sqrt{\frac{1}{2}}\cos(\frac{3\pi}{8})$，则：

<img src="/Users/mac/Library/Application Support/typora-user-images/image-20200621142737327.png" alt="image-20200621142737327" style="zoom:33%;" />

这里我们可以看到$a$，$b$，$c$都是实数，而我们需要转换的像素点是整数。对于实数$DCT$变换，由于在解码端浮点运算的精度问题，会造成解码后数据的失配，引起飘移（见附录1）。所以必须要对标准的$4\times4 DCT$变换改造，采用整数$DCT$变换，减少计算量，同时不损失图像的准确度。可以理解为实数$DCT$变化存在：1. 计算量大。2. 精度问题，导致$H264$中采用的是整数$DCT$变换。

<img src="/Users/mac/Library/Application Support/typora-user-images/image-20200621144100391.png" alt="image-20200621144100391" style="zoom:33%;" />

其中$d=c/b$约等于$0.414$，$\theta$表示的是对位乘。为了简化计算$d=0.5$，同时为了保证正交性，对$b$进行修正，$b=\sqrt\frac{2}{5}$。对矩阵$C$的第二行和第四行，$C^T$的第二列和第四列元素乘以2，相应的改造$E$--->$E_f$，得到：

<img src="/Users/mac/Library/Application Support/typora-user-images/image-20200621145028304.png" alt="image-20200621145028304" style="zoom:33%;" />

对于$\theta$运算，将其归纳到**量化运算**中。所以实际上$DCT$变换只剩下$C_fXC_f^T$，只剩下加法，减法和移位操作。上式运算结果与通常的$DCT$变换结果近似。

对于实际的$DCT$变换$C_fXC_f^T$中的矩阵乘法可以改造为两次的$DCT$变换，首先对行，然后是列。而每一次的整数$DCT$都可以采用蝶形快速算法，节省计算时间。

<img src="/Users/mac/Library/Application Support/typora-user-images/image-20200621145849190.png" alt="image-20200621145849190" style="zoom:33%;" />

## 附录

1. 标准的$DCT$变换是否**无损**可逆？

   关于这个问题，使用了一维$DCT$变换中的$DCT$和$IDCT$公式进行测试，计算的结果如下所示：

   <img src="/Users/mac/Library/Application Support/typora-user-images/image-20200621132529208.png" alt="image-20200621132529208" style="zoom:50%;" />

   从$IDCT$的结果和原始的比较看，由于是浮点运算所以结果和原始的还是有一些出入的，所以说应该是无损的。

2. 二维$DCT$是否可以使用一维$DCT$变换实现？

   通过测试发现二维$DCT$变换和一维$DCT$变换的结果是不一致的。**两者不可以混为一谈**。

3. 正交矩阵

   $AA^T=E$，$E$是单位矩阵，则$A$是正交矩阵。