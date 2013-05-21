## Introduction

This is a very simple 1D colocated shallow water solver. I am just experimenting different methods in this code.

## Build

This program depends on `libjsoncpp`. You can install it on ubuntu by

```
sudo apt-get install libjsoncpp-dev
```

and then use `make` command to compile.

## Governing equations

Governing equations are momentum and continuety equations:

$\frac{\partial p}{\partial t} + u \frac{\partial p}{\partial x} + gh\frac{\partial \xi}{\partial x}=0$ (1)

$\frac{\partial \xi}{\partial t} + \frac{\partial p}{\partial x} =0$ (2)

where $h=d+\xi$, $p=uh$ and $d$ is base depth.

## Numerical method

### Fractional stepping

We follow Kim and Choi (2000):

$\frac{p^*-p^n}{\Delta t} + \frac{F^c_e(U^n,p)-F^c_w(U^n,p)}{\Delta x} =0$ (3)

where $F^c$ is convective flux using FROMM scheme and $U$ is the face velocity.

Remaining parts of momentum equation:

$\frac{p^{n+1}-p^*}{\Delta t} + gh\frac{\partial \xi}{\partial x} =0$ (4)

Taking divergence of above equation and applying continuety on it gives:

$\xi^{n+1} = \xi^n + g (\Delta t)^2 \frac{\partial}{\partial x}\left( h\frac{\partial \xi}{\partial x}\right) - \Delta t \frac{\partial p^*}{\partial x}$ (5)

Overall algorithm is as follows:

1. calculate $p^*$ from Eq. 3
2. calculate $\xi^{n+1}$ from Eq. 5
3. calculate $p^{n+1}$ from Eq. 4
4. update $u$ and $h$

### Momentum interpolation method

To eliminate pressure checkerboarding, $U$ will be computed using:

$U= \left( \frac{p^*}{h} \right) \bar - g \Delta t \frac{\patial \xi}{\partial x}$