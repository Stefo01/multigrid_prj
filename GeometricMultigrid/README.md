# Geometric Multigrid implementation

### Compilation

* From GeometricMultigrid folder, create a build folder build
```
GeometricMultigrid $ mkdir build
```
* Move inside builder folder and then run "make" cmd to compile the entire codebase.
```
GeometricMultigrid $ cd build
build $ cmake ..
build $ make
```
* Navigate inside bin folder to find Multigrid executable and execute it with your custom parameters like:
```
build $ cd bin
bin $ ./Multigrid -n 385 -a 1 -w 10 -ml 5 -test 1 -smt 2
```

### Parallel Compilation

For compiling a parallel version of the program, you have to use, instead of simple cmake cmd:

``` bash
~$ cmake .. -DBUILD_PARALLEL=ON
```

This command will compile the code with specific flags, OpenMP flag, to enable parallel processing. After that one, you have to execute:

```bash
build $ make
build $ cd bin
bin $ export OMP_NUM_THREADS=2
bin $ ./Multigrid -n 385 -a 1 -w 10 -ml 5 -test 1 -smt 2
```

You can set the OMP_NUM_THREADS value on whatever number of parallel threads you want.

### Execution Parameters:

To explore available options and input parameters, utilize the '--help' command. This will provide comprehensive information on the inputs that can be provided to the program. But, in general, you'll find:

When executing ./Multigrid, you can provide the following parameters:
* '-a' (or --alpha):
    Specifies differential constant.
* '-w' (or --width):
    Specifies the width of the squared domain.
* '-n' (or --intervals):
    Sets the discretization number used to discretize the domain.
* '-ml' (or --multigrid):
    Specifies levels of multigrid
* '-test':
    Specifies test for the exe
* '-smt':
    Specifies smoother algorithm, choosing between Jacobi, BiCGSTAB, and Gauss-Seidel


### Test

To conduct a test using an Ubuntu machine, navigate to the "/test/test.ipynb" directory. Modify the desired parameters, including the width and length of the squared domain, discretization number, level of multigrid coarseness and the number of test functions. Choose from the following available test functions:

* Test 0:

$$f(x,y) = 1$$
$$g(x,y) = 0$$

* Test 1:

$$f(x, y) = -5.0 \cdot e^{x} \cdot e^{-2.0 \cdot y}$$
$$g(x, y) = e^{x} \cdot e^{-2.0 \cdot y}$$


* Test 2:

$$ f(x, y) = \sin(k \cdot \sqrt{x^2 + y^2})$$
$$g(x, y) = -k \left( \frac{\cos(kr)}{r} - k\sin(kr) \right)$$


Upon completion, refer to the "result" file to examine the matrix solution. Additionally, a 3D graphic representation of the solution is available for visualization inside the python notebook. 
Note that f is the function defined within the domain, while g specifies the function on the boundary.

# Run with User Interface

---

## Step 1: Verify PHP installation on your localhost

From main folder, run cmd:

```bash
~$ php -v
```

if it reports an error, install php in your linux machine with:

```bash
~$ sudo apt update
~$ sudo apt install php php-cli php-fpm php-xml php-mbstring
```

then re-run to check the installation: 

```bash
~$ php -v
```

After that, navigate to the folder WebInterface and start the local php server:

```bash
GeometricMultigrid $ cd WebInterface
GeometricMultigrid $ php -S localhost:8000
```

Then access the interface via browser at the url: http://localhost:8000/home.php

# Run with python notebook

Navigate to the test folder, where you'll find the 'test.ipynb' file. This file contains all the code you need to run and visualize your results in a clear graphical interface.