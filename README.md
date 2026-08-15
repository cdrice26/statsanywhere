# StatsAnywhere
StatsAnywhere is an experimental, incomplete statistics & scientific computing library implemented in C. Currently, it supports some common linear algebra computations, special functions, and probability distribution functions. Special functions are computed using Gaussian quadrature.

## Why Do This When So Many Scientific Libraries Already Exist?
The intention was mostly for this to be a learning project into the underpinnings of probability distributions as well as low-level coding. It also refreshed some of my linear algebra and calculus knowledge along the way. It was honestly a fun challenge - while I did ask AI to get me unstuck a few times (and check my work), particularly in how to use the Golub-Welsch algorithm and in converting from a general QR decomposition process to a tridiagonal one, by the end of it I was doing most of the work independently and just using it to check my work for silly mistakes, so it's safe to say I learned a fair bit from this project! The only actual practical reason outside personal learning to create this is that implmenting it in C allows it to be imported by virtually any programming language, adding scientific and statistical capabilities to languages that currently lack robust libraries for these tasks.

## Should I Use It?
Probably not. It's not nearly as battle-tested as most scientific libraries, and some functions produce noticeably less precise output. Additionally, since the goal was just to get probability functions working, it's not really a full scientific library anyway. Feel free to download it and play around with it, but it's definitely not production-ready.

## Is It Documented?
Yes, the library has documentation comments and a complete API documentation is available [here](./API.md).

## How Do I Build It?
StatsAnywhere is very simple to build. Simply run the `compile.sh` file and then include `statsanywhere.h` in your project. Be sure to link with the compiled output file.
