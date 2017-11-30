#ifndef PHIST_TPETRA_SDMAT_HPP
#define PHIST_TPETRA_SDMAT_HPP

#include "Kokkos_Core.hpp"
#include "Kokkos_View.hpp"
#include "Tpetra_MultiVector_decl.hpp"
#include "Teuchos_RCP.hpp"

#include <utility>

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
template <typename MV>
class SmallDenseMatrix
{
    typedef typename MV::scalar_type scalar_type;
    typedef typename MV::dual_view_type dual_view_type;
    typedef typename MV::local_ordinal_type local_ordinal_type;
    typedef typename MV::node_type node_type;
    typedef typename MV::device_type device_type;
    typedef typename MV::map_type map_type;

        // raw data
        Teuchos::RCP<dual_view_type> data_;        
        // MV that views the raw data
        Teuchos::RCP<MV> tpetraView_;
        std::pair<int, int> dimension_;

public:
        SmallDenseMatrix(int nrows, int ncols);
        SmallDenseMatrix(Teuchos::RCP<dual_view_type> data,
                         Teuchos::RCP<MV> tpetraView);
        
        // Delete copy/move constructors/assignment operators
        SmallDenseMatrix(SmallDenseMatrix const& sdm) = delete;
        SmallDenseMatrix(SmallDenseMatrix&& sdm) = delete;
        SmallDenseMatrix& operator=(SmallDenseMatrix const& sdm) = delete;
        SmallDenseMatrix& operator=(SmallDenseMatrix&& sdm) = delete;
        
        //void putValue(scalar_type const& scalar);
        auto subview(int rmin, int rmax,
                     int cmin, int cmax);
        auto subcopy(int rmin, int rmax,
                     int cmin, int cmax) const;
        //void update();
        //void multiply();
};

template <typename MV>
SmallDenseMatrix<MV>::SmallDenseMatrix(int nrows, int ncols)
:
    dimension_{std::make_pair(nrows, ncols)}
{
    // allocate nrows x ncols matrix
    auto dualView = Kokkos::View<scalar_type**>("sdMat", nrows, ncols);
    data_ = Teuchos::rcp(dualView, false);
}

template <typename MV>
SmallDenseMatrix<MV>::SmallDenseMatrix(Teuchos::RCP<dual_view_type> data,
                                       Teuchos::RCP<MV> tpetraView)
:
    data_{data},
    tpetraView_{tpetraView}
{}

template <typename MV>
auto SmallDenseMatrix<MV>::subview(int rmin, int rmax,
                                   int cmin, int cmax)
{
    using std::make_pair;
    // phist wants inclusive endpoints, hence the + 1
    auto view = Kokkos::subview(data_, make_pair(rmin, rmax + 1), 
                                       make_pair(cmin, cmax + 1));
    return Teuchos::rcp(new MV(tpetraView_->getMap(), view), false);
}

template <typename MV>
auto SmallDenseMatrix<MV>::subcopy(int rmin, int rmax,
                               int cmin, int cmax) const
{
    using std::make_pair;
    // phist wants inclusive endpoints, hence the + 1
    const auto view = Kokkos::subview(data_, make_pair(rmin, rmax + 1), 
                                             make_pair(cmin, cmax + 1));
    auto copy = Teuchos::rcp(new MV(tpetraView_->getMap(), view), false);
    // Deep copy the data into the new MV
    copy->assign(*view);
    return copy;
}


}// namespace tpetra
}// namespace phist
#endif
