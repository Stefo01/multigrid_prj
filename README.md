# Implementation of Multigrid method

This project presents a comprehensive study and implementation of multigrid methods, as part of the Exam for of the Advanced Methods for
Scientific Computing course at Politecnico di Milano.

## Method description

The Multigrid method is an iterative numerical technique used to solve discretized systems of equations, particularly those arising from partial differential equations (PDEs). It's a powerful method for accelerating the convergence of iterative solvers, especially for large-scale problems.

Multigrid methods can be broadly categorized into two principal classes:

* Geometric Multigrid (GMG)
* Algebraic Multigrid (AMG)

You'll find the implementation of Geometric Multigrid within the `GeometricMultigrid` folder in this repository, while Algebraic Multigrid can be found and used in the `AMG` folder. 

We've chosen to keep these two methods separate due to differences in their compilation processes and initialization, allowing us to better focus on each method's unique capabilities. This separation also highlights a fundamental distinction in how they handle the problem domain:

* The **Geometric Multigrid** (GMG) implementation is designed to build its own mesh and then execute the program, handling the entire discretization process internally.
* In contrast, the **Algebraic Multigrid** (AMG) implementation focuses its code on operating with an existing `.msh` mesh file, which you can place inside the mesh folder. This means AMG is more adaptable to pre-existing mesh structures.

You can find the instruction to run and execute AMG and GMG inside their readme: 

* [AMG](./AMG/README.md)
* [GMG](./GeometricMultigrid/README.md)

## License

[MIT License](LICENSE) (or specify your license)

---

**Contact:**  
For questions or contributions, please open an issue or pull request.