# Trying to get a glimpse what we have to do

- navier stokes flow around ball in a tube (this is constant, only one solve used)
- shape derivative to change shape of ball
  - cost function: drag (minimize)
- iterate over this
- problems:
  - make stokes fem easily available with fitting functions for shape derivative
  - shape derivative maths overall...?
  - discretization between shape of ball in mesh and some function for shape derivative?
    - is this actually even needed?
- initial change function? We have to start with something to change on the grid?!?!
- Change functions on our grid must have compact support - boundaries can't change, only ball (tube is static)

## things to look into/I should learn first
- state equation
- regression of 2d function
- cost functions for shape derivatives

## What should actually be happening
##### understood a step? put it in here
- navier stokes flow around ball on fixed domain
-

## available tutorials/examples we have to fix together + their explanation

#### 1. [navier stokes flow example](https://docu.ngsolve.org/latest/i-tutorials/unit-3.2-navierstokes/navierstokes.html) (def used and basically done)

- unsteady navierstokes equation
- "Taylor-Hood velocity-pressure pair"
  - we need both velocity and pressure field to calculate this
- stokes equation ngsolve: `stokes = (nu*InnerProduct(grad(u),grad(v))-div(u)*q-div(v)*p)*dx`
- don't need the time discetization -> we just need one solve to get current drag
``` python
drag_x_test = GridFunction(X)
# test function for velocity on opposite direction
drag_x_test.components[0].Set(CoefficientFunction((-20.0,0)), definedon=mesh.Boundaries("cyl"))
# restoring initial data
res.data = f.vec - a.mat*gfu.vec
drag = InnerProduct(res, drag_x_test.vec)
```

### different shape derivative thingies
#### 2. [Example](https://docu.ngsolve.org/latest/i-tutorials/unit-7-optimization/01_Shape_Derivative_Levelset.html) with previous given function that is negative at goal points (prob useless)
- given function that is negative where the final shape should exist (and positive where it shouldn't)
- start with normal circle as shape and go for iterations
- somehow they end up with this formula for the shape derivative problem:
<br><img src="https://render.githubusercontent.com/render/math?math=\LARGE\color{white}DJ(\Omega)(X) = \int_\Omega f {div(X)} - - \nabla f\cdot X\dx">
> plus doesnt work in equation??

##### single optimization step:
- calculate with stepsize alpha:
<br><img src="https://render.githubusercontent.com/render/math?math=\LARGE\color{white}\Omega_1 = ({id} - \alpha_1 X_0)(\Omega_0)">
- test one deformation, calculate next change out of it
  - one deformation for aX and dJOmega (see below, used for first calculation afterwards)
``` python
mesh.SetDeformation(gfset) # deform current domain with gfset
aX.Assemble() # assemble on deformed domain
dJOmega.Assemble()
mesh.UnsetDeformation() # unset deformation
```


- X0 (initial guess) calculated by `gfX.vec.data = aX.mat.Inverse(bunch of parameters) * dJOmega.vec`
  - aX: Bilinear form: (where does this come from, why is it chosen like this?)
<br><img src="https://render.githubusercontent.com/render/math?math=\Large\color{white}(\varphi,\psi) \mapsto b_\Omega(\varphi,\psi):= \int_\Omega (\nabla \varphi+\nabla \varphi^\top): \nabla\psi -- \varphi \cdot \psi\dx">
  - dJOmega: Shape derivatve (see first equation in this section)
  - code for this:
  - f in here is the funny function that is negative where we want the shape
``` python
VEC = H1(mesh, order=1, dim=2) # vector space for shape gradient
gfset = GridFunction(VEC) # grid function for deformation field
gfX = GridFunction(VEC)
PHI, PSI = VEC.TnT() # Test and trial functions
dJOmega = LinearForm(VEC) # shape derivative
dJOmega += (div(PSI)*f + InnerProduct(grad_f, PSI) )*dx
```
- cost: `Integrate(f*dx, mesh)`
- stepsize: `alpha = 20.0 / Norm(gfX.vec)`
- changing the mesh accordingly:
  - gfset changed with the calculated shape derivative
  - gfset used on mesh -> cost goes down (hopefully)
``` python
gfset.vec[:]=0
gfset.vec.data -= alpha * gfX.vec
mesh.SetDeformation(gfset)#
new_cost = Integrate(f, mesh)
mesh.UnsetDeformation()
```
- mesh is always set back to its original, only gfset is changed/worked on
- faster with L-BFGS method from scipy
- possibility to improve mesh quality via Cauchy-Riemann equations

#### 3. [PDE-Constrained Shape Optimization](https://docu.ngsolve.org/latest/i-tutorials/unit-7-optimization/02_Shape_Derivative_Laplace.html)
- this seems somewhat random/useless
- difference between PDE constrained and the one before?
  - state equation
  - f and ud exist
  - adjoint equation?
- very long/weird shape derivative formula
- rest is pretty similar to example 2

#### 4. [PDE-Constrained Shape Optimization (semi-automated)](https://docu.ngsolve.org/latest/i-tutorials/unit-7-optimization/03_Shape_Derivative_Laplace_SemiAuto.html)
- reference to the paper given to us
- again with state&adjoint equation (whatever that is)
- define Lagrangian, solve with this
- it is possible to write the pde only once and do the entire optimization thing on their own

#### 5. [PDE-Constrained Shape Optimization (fully automated)](https://docu.ngsolve.org/latest/i-tutorials/unit-7-optimization/03b_Shape_Derivative_Laplace_FullyAuto.html)
- step of defining the perturbed Lagrangian can be automated
- solve thing with lagrangian with newton (wasn't done before?!)
- afterwards seems the same as before
