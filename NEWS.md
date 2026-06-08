# cforecast 0.1.1

* New KFAS Backend

Added support for the **KFAS** backend for state-space filtering and smoothing routines. The backend can be selected using the `package` argument (`package = "KFAS"`).

This backend is particularly useful when the forecast error variance matrix is **singular or near-singular**, a situation in which the default **FKF** implementation may fail. While **KFAS** provides a more robust alternative, it is computationally more demanding and may result in longer execution times.
