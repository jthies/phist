//! this functionality allows to create an operator that gives a        
//! different effect for each column of a multi-vector (mvec).          
//! If Y=Op*X, then Y_i = (A + shift_i I) X_i, for an arbitrary         
//! operator A_op. The shifts to be used are set at construction        
//! time, but the order in which they are applied can be changed        
//! by using the functions set/reset_order.                             
//! When the operator is applied to a multi-vector with k columns,      
//! A - sigma[shift_order[i]]*I is applied to column i. If k is larger  
//! than num_shifts, shift 0 is applied to the last k-num_shifts cols.  
typedef struct TYPE(multishift) {
 TYPE(const_op_ptr) A_op_; // data structure needed for computing A*x
 _ST_* shifts_; // stores all possible shifts
 int num_shifts_; // maximum number of shifts the operator can hold
 int* shift_order_;
} TYPE(multishift);



//! create an operator that behaves like a multi-shift operator for a given
//! base operator A_op. The shifts are copied and their given order defines
//! the default order in which they are applied.
void SUBR(msOp_create)(TYPE(op_ptr)* msOp, TYPE(const_op_ptr) A_op, 
        int num_shifts, _ST_* shifts,int* ierr);

//! set entry in the ordering array such that given shift j is applied to vector column i
void SUBR(msOp_set_order)(TYPE(op_ptr) msOp, int i, int j, int* ierr); 

//! reset the order in which the shifts are applied to the default, 0:num_shifts-1
void SUBR(msOp_reset_order)(TYPE(op_ptr) msOp, int* ierr);

//! set entry in the ordering array such that given shift j is applied to vector column i
void SUBR(msOp_delete)(TYPE(op_ptr) msOp, int* ierr); 

