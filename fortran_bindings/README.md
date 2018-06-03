# phist_fort - generate simple Fortran 2003 bindings for [PHIST](https://bitbucket.org/essex/phist)

The phist software can be found at https://bitbucket.org/essex/phist

Based on a PHIST installation, phist_fort runs a preprocessor over the installed headers and generates Fortran modules 
using the iso_c_bindings feature of Fortran 2003. This software is not included in PHIST itself because it is based on
a Python script from the GTK+ library, which is distributed under a GPL license whereas PHIST is distributed under a BSD 
license. The original cfwrapper.py file can be found [here](https://github.com/jerryd/gtk-fortran/blob/master/src/cfwrapper.py)

**DISCLAIMER:** the Fortran bindings are at this stage **largely untested**, we do not make any claims about their 
correctness.

## USAGE EXAMPLE

export CMAKE_PREFIX_PATH=${PHIST_DIR}:${CMAKE_PREFIX_PATH}  
export FC=mpif90  
mkdir build  
cd build  
cmake ..  
make  
make test  

note that you have to tell cmake to use the MPI Fortran compiler, and that 
you have to make it find your phist installation. The final command 'make test',
currently fails because of bugs in the generated interfaces.

## Translation rules applied

In order to create a pretty Fortran interface we need the C code to adhere to some notation conventions.
The data types that appear must be declared in the dictionary (phist_dict.py), along with the f_type they
should yield and the iso_c_binding types that need to be USEd.

Example: c_type="int", f_type="integer(c_int)", iso_c="c_int"

C structs should be written as:  

typedef struct c_struct {  
...  
} c_struct;  

and are treated slightly differently: a pointer to struct is passed as TYPE(c_ptr) in order
to allow passing in C_NULL_PTR (see table below and issue #1)

### struct data members

-------------------------------------------------
|       **C**   |   **Fortran**                 |
|---------------|-------------------------------|
| c_type a;     | f_type :: a                   |
| c_type a[5];  | f_type, dimension(5) :: a     |
-------------------------------------------------

### function arguments

-------------------------------------------------------------------------
|       **C**           |   **Fortran**                                 |
|-----------------------|-----------------------------------------------|
| c_type a              | f_type :: a, value                    :: a    |
| c_type *a             | f_type, intent(inout)                 :: a    |
| const c_type *a       | f_type, intent(in)                    :: a    |
| c_type a[]            | f_type, dimension(*), intent(inout)   :: a    |
| const c_type a[]      | f_type, dimension(*), intent(in)      :: a    |
------------------------|-----------------------------------------------|
| c_struct s            | type(c_struct), value                 :: s    |
| c_struct *s           | type(c_ptr), intent(inout)            :: s    |
| const c_struct *s     | type(c_ptr), intent(in)               :: s    |
| c_struct s[]          | type(c_struct), dimension(*), intent(inout)   |
| const c_struct s[]    | type(c_struct), dimension(*), intent(in)      |
-------------------------------------------------------------------------

In these tables, the order of 'const' and 'c_type' may be exchanged, and the intent(inout) may be omitted by the 
generator.
