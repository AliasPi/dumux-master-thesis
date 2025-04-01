#ifndef SHIFTED_GRID_GEOMETRY_HH
#define SHIFTED_GRID_GEOMETRY_HH

#include <dune/common/fvector.hh>
#include <memory>

// Template wrapper for any BaseGridGeometry.
// It shifts all vertical coordinates (index dim-1) by a constant offset.
template <typename BaseGridGeometry>
class ShiftedGridGeometry : public BaseGridGeometry
{
public:
    using Base = BaseGridGeometry;
    // Instead of Base::Scalar (which does not exist), use the grid view's coordinate type.
    using Scalar = typename Base::GridView::ctype;
    static constexpr int dim = Base::GridView::dimensionworld;
    using FieldVector = Dune::FieldVector<Scalar, dim>;

    // Constructor: forward the grid view to the base class and store the vertical shift.
    ShiftedGridGeometry(const typename Base::GridView& gridView, Scalar verticalShift = 2000)
      : Base(gridView), verticalShift_(verticalShift)
    {}

    // Override the bounding box functions to add the vertical shift.
    FieldVector bBoxMin() const
    {
        FieldVector min = Base::bBoxMin();
        min[dim - 1] += verticalShift_;
        return min;
    }

    FieldVector bBoxMax() const
    {
        FieldVector max = Base::bBoxMax();
        max[dim - 1] += verticalShift_;
        return max;
    }

    // --- Local view wrapper ---
    // The simulation uses local view objects to query subcontrol volumes (SCVs).
    // We wrap these so that the returned dof positions are shifted.

    // Wrapper for a subcontrol volume.
    class ShiftedSCV
    {
    public:
        using SubControlVolume = typename Base::LocalView::SubControlVolume;
        using FieldVector = Dune::FieldVector<Scalar, dim>;

        ShiftedSCV(const SubControlVolume& baseScv, Scalar verticalShift)
          : baseScv_(baseScv), verticalShift_(verticalShift)
        {}

        // Return the dof position shifted by verticalShift.
        FieldVector dofPosition() const
        {
            FieldVector pos = baseScv_.dofPosition();
            pos[dim - 1] += verticalShift_;
            return pos;
        }

        // Provide access to the underlying SCV if needed.
        const SubControlVolume& base() const { return baseScv_; }
    private:
        const typename Base::LocalView::SubControlVolume& baseScv_;
        Scalar verticalShift_;
    };

    // Iterator wrapper to iterate over shifted SCVs.
    template<typename BaseIterator>
    class ShiftedIterator
    {
    public:
        using difference_type   = typename BaseIterator::difference_type;
        using value_type        = ShiftedSCV;
        using pointer           = void;
        using reference         = value_type;
        using iterator_category = std::forward_iterator_tag;

        ShiftedIterator(BaseIterator baseIt, Scalar verticalShift)
          : baseIt_(baseIt), verticalShift_(verticalShift)
        {}

        ShiftedIterator& operator++() { ++baseIt_; return *this; }
        ShiftedIterator operator++(int) { ShiftedIterator tmp(*this); ++baseIt_; return tmp; }
        bool operator==(const ShiftedIterator& other) const { return baseIt_ == other.baseIt_; }
        bool operator!=(const ShiftedIterator& other) const { return baseIt_ != other.baseIt_; }
        value_type operator*() const { return ShiftedSCV(*baseIt_, verticalShift_); }
    private:
        BaseIterator baseIt_;
        Scalar verticalShift_;
    };

    // A wrapper for the local view itself.
    class ShiftedLocalView
    {
    public:
        using LocalView = typename Base::LocalView;
        using Iterator = ShiftedIterator<typename LocalView::Iterator>;

        ShiftedLocalView(const LocalView& baseLv, Scalar verticalShift)
          : baseLv_(baseLv), verticalShift_(verticalShift)
        {}

        // Bind an element (same as in the base local view).
        void bindElement(const typename LocalView::Element& element)
        {
            baseLv_.bindElement(element);
        }

        // Return wrapped iterators over the subcontrol volumes.
        Iterator begin() const { return Iterator(baseLv_.begin(), verticalShift_); }
        Iterator end()   const { return Iterator(baseLv_.end(), verticalShift_); }
    private:
        LocalView baseLv_;
        Scalar verticalShift_;
    };
    // --- End local view wrapper ---

    // Override localView() to return our shifted local view.
    ShiftedLocalView localView() const
    {
        return ShiftedLocalView(Base::localView(), verticalShift_);
    }

private:
    Scalar verticalShift_;
};

#endif // SHIFTED_GRID_GEOMETRY_HH
