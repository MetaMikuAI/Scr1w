# 数论

## 费马小定理

#### 表述

设$a$为整数，$p$为质数，则有：
$$
a^{p} \equiv a \space (mod \space p)
$$
或写作：
$$
a^{p-1} \equiv 1 \space(mod \space p)
$$

#### 证明

##### 证明1（数归）

使用**数学归纳法**进行证明：

- 当$a=1$时，显然有

$$
a^{p} = a(mod \space p)
$$

- 设$a=k$时原式成立，则$a=k+1$时，有：

$$
(a+1)^{p} &\equiv& \sum_{k=0}^{p}\binom{p}{k}a^{k} &\space(mod \space p)&\\
&\equiv& \sum_{k=0}^{p}\frac{p!}{k!(p-k)!} &\space(mod \space p)&\\
&\equiv& a+0+ \cdots +0+1 &\space(mod \space p)&\\
&\equiv& a+1 &\space(mod \space p)&
$$

即$a=k+1$仍成立

- 综上，对任意$a\in \Z^+$，有：$a^{p} = a(mod \space p)$成立。

证毕。

### 相关资料：

1. [wiki](https://zh.wikipedia.org/zh-cn/%E8%B4%B9%E9%A9%AC%E5%B0%8F%E5%AE%9A%E7%90%86)

## 欧拉$$\varphi$$函数定理

~~突然发现了可以行间插入latex，马上去改~~

### 表述

如果$n,a \in \Z^+$且$gcd(a,n)=1$，则有：
$$
a^{\phi(n)} \equiv 1 \space (mod \space n)
$$

### 证明

好像要用到群论，先挂上

## RSA

### 流程

1. 取两个大质数$$p,q$$
2. 取$$n=p \times q$$，$$\phi=(p-1)\times (q-1)$$
3. 取合适的私钥$$d$$,公钥$$e$$使$$d\times e \space mod \space \phi = 1$$
4. 销毁$$p,q$$
5. 公开模数$$n$$,公钥$$e$$
6. 加密$$c=m^{e} \space mod \space n$$
7. 传回密文$$c$$
8. 解密$$m=c^{d} \space mod \space n$$

###  解密原理

$$
m &\equiv& c^{d} &\space (mod \space n)\\
&\equiv& (m^{e} \space mod \space n)^{d} &\space (mod \space n)\\
&\equiv& (m^{e})^{d} &\space (mod \space n)\\
&\equiv& m^{de} &\space (mod \space n)\\
&\equiv& m^{de \space mod \space \phi} &\space (mod \space n)\\
&\equiv& m &\space (mod \space n)\\
$$

###  安全性

## ECC

直接去看方法，真的没看懂，有很多问题：

- 圆锥曲线计算难道不会出现超级多的浮点数吗？在大量运算之下可以保证精度吗？
- 如果用模数的话，怎么能模出合理的数？或者模后点变离散后，要是A+B=不存在的点呢？
- 无穷远点是哪条直线的极点？？难道不是无穷个无穷远点吗？？哪条符合规定

所以决定直接去看代码实现

在知乎上找到了一篇文章[ECC椭圆曲线密码学的原理、公式推导、例子、Python实现和应用](https://zhuanlan.zhihu.com/p/42629724)

```python
# ECC在Fp域上的点集Python实现
def show_points(p, a, b):
    return [(x, y) for x in range(p) for y in range(p) if (y*y - (x*x*x + a*x + b)) % p == 0]

print(show_points(p=23, a=1, b=1))
#[(0, 1), (0, 22), (1, 7), (1, 16), (3, 10), (3, 13), (4, 0), (5, 4), (5, 19), (6, 4), (6, 19), (7, 11), (7, 12), (9, 7), (9, 16), (11, 3), (11, 20), (12, 4), (12, 19), (13, 7), (13, 16), (17, 3), (17, 20), (18, 3), (18, 20), (19, 5), (19, 18)]
```

啊？直接把$\R^2$拆成了$(\Z\bigcap\left [ 0,p \right))^2$？？或者说是$E \bigcap (\Z\bigcap\left [ 0,p \right))^2$?~~（这里后来发现错了）~~

下一段代码：

```python
# ECC在Fp域上加法、倍乘运算

# 求value在Fp域的逆——用于分数求逆
def get_inverse(value, p):
    for i in range(1, p):
        if (i * value) % p == 1:
            return i
    return -1

# 求最大公约数——用于约分化简
def get_gcd(x, y):
    if y == 0:
        return x
    else:
        return get_gcd(y, x % y)  

# 计算P+Q函数
def calculate_p_q(x1,y1,x2,y2,a,b,p):
    flag = 1  # 控制符号位
    
    # 若P = Q，则k=[(3x1^2+a)/2y1]mod p
    if x1 == x2 and y1 == y2:
        member = 3 * (x1 ** 2) + a  # 计算分子
        denominator = 2 * y1        # 计算分母

    # 若P≠Q，则k=(y2-y1)/(x2-x1) mod p
    else:
        member = y2 - y1
        denominator = x2 - x1 
        if member* denominator < 0:
            flag = 0
            member = abs(member)
            denominator = abs(denominator)
    
    # 将分子和分母化为最简
    gcd_value = get_gcd(member, denominator)
    member = member // gcd_value
    denominator = denominator // gcd_value

    # 求分母的逆元    
    inverse_value = get_inverse(denominator, p)
    k = (member * inverse_value)
    if flag == 0:
        k = -k
    k = k % p

    # 计算x3,y3
    """
        x3≡k^2-x1-x2(mod p)
        y3≡k(x1-x3)-y1(mod p)
    """
    x3 = (k ** 2 - x1 - x2) % p
    y3 = (k * (x1 - x3) - y1) % p
    return [x3,y3]
    
# 计算nP函数
def calculate_np(p_x, p_y,a,b,p):
    tem_x = p_x
    tem_y = p_y
    p_value = calculate_p_q(tem_x,tem_y, p_x, p_y,a,b,p)
    tem_x = p_value[0]
    tem_y = p_value[1]
    return p_value
        
# PPT例2、例3：y^2=x^3+x+1(mod 23)
p, a, b = 23, 1, 1

# PT例2计算P+Q ,其中p=(3,10)、q=(9,7)
print(calculate_p_q(3,10,9,7,1,1,23))
 
 # PT例3计算2P ,其中p=(3,10)
print(calculate_np(3,10,1,1,23))
```

里面有一些函数，根据前的学习，我改写一下：

```python
# ECC在Fp域上加法、倍乘运算
from gmpy2 import gcd
# 求value在Fp域的逆——用于分数求逆
def get_inverse(value, p):
    test = pow(value,-1,p)
    for i in range(1, p):
        if (i * value) % p == 1:
            return i
    return -1

# 计算P+Q函数
def calculate_p_q(x1,y1,x2,y2,a,b,p):
    flag = 1  # 控制符号位
    
    # 若P = Q，则k=[(3x1^2+a)/2y1]mod p
    if x1 == x2 and y1 == y2:
        member = 3 * (x1 ** 2) + a  # 计算分子
        denominator = 2 * y1        # 计算分母

    # 若P≠Q，则k=(y2-y1)/(x2-x1) mod p
    else:
        member = y2 - y1
        denominator = x2 - x1 
        if member* denominator < 0:
            flag = 0
            member = abs(member)
            denominator = abs(denominator)
    
    # 将分子和分母化为最简
    gcd_value = gcd(member, denominator)
    member = member // gcd_value
    denominator = denominator // gcd_value

    # 求分母的逆元    
    inverse_value = get_inverse(denominator, p)
    k = (member * inverse_value)
    if flag == 0:
        k = -k
    k = k % p

    # 计算x3,y3
    """
        x3≡k^2-x1-x2(mod p)
        y3≡k(x1-x3)-y1(mod p)
    """
    x3 = (k ** 2 - x1 - x2) % p
    y3 = (k * (x1 - x3) - y1) % p
    return [x3,y3]
    
# 计算nP函数
def calculate_np(p_x, p_y,a,b,p):
    tem_x = p_x
    tem_y = p_y
    p_value = calculate_p_q(tem_x,tem_y, p_x, p_y,a,b,p)
    tem_x = p_value[0]
    tem_y = p_value[1]
    return p_value
        
# PPT例2、例3：y^2=x^3+x+1(mod 23)
p, a, b = 23, 1, 1

# PT例2计算P+Q ,其中p=(3,10)、q=(9,7)
print(calculate_p_q(3,10,9,7,1,1,23))
 
 # PT例3计算2P ,其中p=(3,10)
print(calculate_np(3,10,1,1,23))
```

改写的时候了解了大致逻辑，用matlab跑图看看

以$y^2=x^3+x+1,p=23$为例：

```matlab
ezplot('x^3+1*x+1-y^2',[-23,23,-23,23])
```



![ECC1](E:\Desktop\MetaMiku\Scr1w\SSSCTF\images\ECC1.png)

然而发现例子中的点$(3,10)$并不在这条线上，重新回到show_points中仔细看了下，发现对点的判定为：

```python
if (y*y - (x*x*x + a*x + b)) % p == 0
#                            ^^^
```

啊？

继续绘图

```matlab
for a = -20:1:20
	ezplot(['x^3+1*x+1-y^2+23*',num2str(a)],[-23,23,-23,23])
	hold on
end
```

![ECC2](E:\Desktop\MetaMiku\Scr1w\SSSCTF\images\ECC2.png)

乱死了，不取模的话好多非必要的点，重画

```matlab
p=23;
a=1;
b=1;
for x = 0:1:p-1
	for y = 0:1:p-1
		if mod(y^2 - (x^3 + a*x + b),p)==0
			plot(x,y,'b*')
			hold on
		end
	end
end
```

![ECC3](E:\Desktop\MetaMiku\Scr1w\SSSCTF\images\ECC3.png)

好看多了，但是看不出来是圆锥曲线了（

图像关于$x$轴对称

~~我要是会manim就好了~~

以$P=Q=(3,10)$为例子，计算$P+Q$中的$k$：（$P \ne Q$时的$k$要更好理解，故略）

先计算切线
$$
y^2 = x^3 + ax + b\\
2y \mathrm{d}y = 3x^2 \mathrm{d}x+a \mathrm{d}x\\
k = \frac{\mathrm{d}y}{\mathrm{d}x}=\frac{3x^2+a}{2y}\\
k \equiv \frac{3x^2+a}{2y} \pmod{p}
$$
当$a=1,x=3,y=10,p=23$时，$k=6$

接下来用$k$求$P+Q$，暂时忽略$\pmod{p}$
$$
\begin{cases}
PQ: y=k(x-x_1)+y_1\\
E: y^2 = x^3 + ax + b
\end{cases}
$$
将$PQ$带入$E$：
$$
k^2(x-x_1)^2+y_1^2+2ky_1(x-x_1) = x^3 + ax + b\\
k^2x^2-2k^2x_1x+k^2x_1^2+y_1^2+2ky_1x-2ky_1x_1 = x^3 + ax + b\\
x^3 -k^2x^2 +(2k^2x_1-2ky_1+a)x+ (b -k^2x_1^2 -y_1^2+-2ky_1x_1)=0\\
$$
待定系数法：
$$
(x-x1)(x-x2)(x-x3)=0\\
x^3-(x_1+x_2+x_3)x^2+(x_1x_2-x_1x_3-x_2x_3)x+x_1x_2x_3=0
$$
对应项系数相等：
$$
\begin{cases}
x_1+x_2+x_3 =  k^2\\
2k^2x_1-2ky_1+a=2k^2x_1-2ky_1+a\\
b -k^2x_1^2 -y_1^2+-2ky_1x_1 = x_1x_2x_3
\end{cases}
$$

易得$x_3$：
$$
x_3 = k^2-x_1+x_2
$$
带入直线$PQ$：
$$
y_3=k(x_3-x_1)+y_1
$$
由于需要取关于$x$轴的对称点，实际$y_3$应为：
$$
y_3 = k(x_1-x_3)-y_1
$$
综上：
$$
\begin{cases}
x_3 = k^2-x_1+x_2\\
y_3 = k(x_1-x_3)-y_1
\end{cases}
$$
重构ECC示例脚本：

```PYTHON
import gmpy2
def sgn(x):
    if x > 0:
        return 1
    if x < 0:
        return -1
    if x == 0:
        return 0
class frac:
    def __init__(self,sign,numerator,denominator):
        self.sign = sign  
        self.numerator = numerator
        self.denominator = denominator
    def simplify(self):
        temp = frac(self.sign,self.numerator,self.denominator)
        if sgn(temp.denominator)*sgn(temp.numerator)==-1:
            temp.sign=-temp.sign
        temp.numerator = abs(temp.numerator)
        temp.denominator = abs(temp.denominator)
        gcd_value = gmpy2.gcd(temp.numerator, temp.denominator)
        return frac(temp.sign,temp.numerator // gcd_value , temp.denominator // gcd_value)
    def invert(self,p):
        return (self.sign*self.numerator * pow(self.denominator,-1,p)) % p
    def print(self):
        print(('+ 'if self.sign==1 else '- ')+str(self.numerator)+' / '+str(self.denominator))
class Point:
    def __init__(self,x,y):
        self.x = x
        self.y = y
    def __eq__(self,other):
        if self.x == other.x and self.y == other.y:
            return True
        else:
            return False
    def __add__(self,other):
        global E
        if self==other:
            k=frac(1, 3*(self.x**2)+E.a , 2*self.y)
        else:
            k = frac(1, other.y-self.y , other.x-self.x)
        k = k.simplify()
        k = k.invert(E.p)
        x3 = (k ** 2 - self.x - other.x) % E.p
        y3 = (k * (self.x - x3) - self.y) % E.p
        R = Point(x3,y3)
        return R
    def __mul__(self,n):  #只能写成右乘数了T_T
        T = Point(self.x,self.y)
        R = Point(self.x,self.y)
        n -= 1
        while n:
            if n & 1:
                R = R+T
            T = T+T
            n = n>>1
        return R
    def print(self):
        print((self.x,self.y))  
class Ellipse:
    def __init__(self,a,b,p):
        self.a=a
        self.b=b
        self.p=p

E = Ellipse(1,1,23)
P = Point(3,10)
(P+P).print()
(P*114514).print()
Q = Point(13,7)
(P+Q).print()
(Q+P).print()
```

##### 回答一下开始的疑问

- 直观的讲，如果没有模运算的话，这些都是有理数，可以写成分数，而有了分数就可以对分母取模运算逆元，这样将有理数域压缩到了整数域

- 算出不存在的点，那就是无穷远点啦，考虑两种情况：

- 1. $PQ$斜率存在，则对于椭圆曲线某一点的切线斜率：

  $$
  y^2 = x^3 + ax + b\\
  2y \mathrm{d}y = 3x^2 \mathrm{d}x+a \mathrm{d}x\\
  k = \frac{\mathrm{d}y}{\mathrm{d}x}=\frac{3x^2+a}{2y}\\
  k^2 = \frac{9x^4+a^2+6ax^2}{4y^2}\\
  k^2 = \frac{9x^4+a^2+6ax^2}{4x^3+4ax+4b}\\
  \lim_{x \to \infty }{k^2} = \lim_{x \to \infty }{\frac{9x^4+a^2+6ax^2}{4x^3+4ax+4b}}=\infty
  $$

  可见，椭圆曲线远离$y$轴的部分斜率趋于正负无穷，也就是只要$PQ$斜率存在，总能和椭圆曲线$E$交上

  2. $PQ$斜率不存在，即与$x$轴垂直，与$E$有且仅有两个交点(或重合的两个交点)，要有第三个交点，那就是无穷远点$O(0,1,0)$*~~（喜闻乐见的齐次坐标）~~*

## 辗转相除法/欧几里得算法/Euclidean algorithm

##### 核心原理

> 两个数的最大公约数等于其中较小的数字和二者之间余数的最大公约数

##### 证明（非严谨证明）

若

$$
M=ax,N=bx(M>N)\\
gcd(a,b)=1,gcd(M,N)=x
$$

则
$$
M \space mod \space N = (a-[\frac{M}{N}]b)x\\
gcd(a-[\frac{M}{N}]b,b)=1,gcd(M \space mod \space N,N) = x
$$

##### 代码示例

```C
int gcd(int m, int n){
    int r = m%n;
    while (r!=0){
		m = n;
		n = r;
        r = a%b;
    }
    return n;
}
```

##### 相关资料

[漫画算法：辗转相除法是什么鬼？](https://zhuanlan.zhihu.com/p/31824895)

