link_libraries(mps gmp m)

add_executable(mandelbrot-solver
  mandelbrot-poly.h
  mandelbrot-evaluation.c 
  mandelbrot-poly.c 
  mandelbrot-newton.c 
  mandelbrot-solver.c 
  mandelbrot-starting.c)
