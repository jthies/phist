#ifndef PHIST_TPETRA_SDMAT_HPP
#define PHIST_TPETRA_SDMAT_HPP

namespace phist {
namespace tpetra {

//! class for representing a small and dense (node local) matrix in Tpetra/Kokkos

//! This class is a light-weight version of Tpetra::MultiVector. It's central data object
//! is a Kokkos DualView with a host and a device side (if applicable). As soon as an arithmetic
//! function (a member of Tpetra::MultiVector) is needed, we create a Tpetra Map and wrap the
//! data as a full-blown Tpetra::MultiVector. Hence, just using the sdMat for views and copies,
//! setting values etc. does not require the map to be created and should be cheap.
//!
//! this class is templated on the MultiVector type it should be compatible with.
template<class MV>
class SmallDenseMatrix
{

    typedef MV::scalar_type scalar_type;
    typedef MV::dual_view_type dual_view_type;
    typedef MV::local_ordinal_type local_ordinal_type;
    typedef MV::node_type node_type;
    typedef MV::device_type device_type;
    typedef MV::map_type map_type;


// data members
private:

        Teuchos::RCP<dual_view_type> data_;
        
        Teuchos::RCP<MV> tpetraView_;
};
}// namespace tpetra
}// namespace phist
#endif
