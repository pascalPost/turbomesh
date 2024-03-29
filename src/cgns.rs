// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

pub mod interface {
    use cgns_sys::*;
    use std::ffi::{CStr, CString};

    pub use std::error::Error;
    pub use std::os::raw::{c_int, c_void};

    /// rust wrapper for CG_MODE_* (excluding CG_MODE_CLOSED which is not excepted
    /// for the cg_open function)
    pub enum CgMode {
        _Read,
        Write,
        _Modify,
    }

    impl CgMode {
        /// return the CG_MODE_* constant
        fn value(&self) -> u32 {
            match self {
                CgMode::_Read => CG_MODE_READ,
                CgMode::Write => CG_MODE_WRITE,
                CgMode::_Modify => CG_MODE_MODIFY,
            }
        }
    }

    /// checks the cgns function's returned integer error codes and wraps the
    /// error massage into the Result
    ///
    /// From the cgns doc:
    /// All C functions return an integer value representing the error status.
    /// An error status different from zero implies that an error occured.
    /// The error message can be printed using the error handling functions of the
    /// CGNS library.
    /// The error codes are coded in the C and Fortran include files cgnslib.h and
    /// cgnslib_f.h.
    fn check(ier: c_int) -> Result<(), Box<dyn Error>> {
        if ier == 0 {
            Ok(())
        } else {
            unsafe {
                let err = CStr::from_ptr(cg_get_error()).to_str()?;
                Err(format!("cgns error encountered : {}.", err))?
            }
        }
    }

    /// cgns structure name and labels char length checker
    fn check_str_len(s: &str) -> Result<(), &str> {
        if s.len() > 32 {
            Err("CGNS structure names and labels limited to only 32 characters.")?
        }

        Ok(())
    }

    /// Open a CGNS file.
    ///
    /// Arguments:
    /// * `filename` : Name of the CGNS file, including path name if necessary. There is no limit on the length of this character variable.
    /// * `mode` : Mode used for opening the file. The modes currently supported are CG_MODE_READ, CG_MODE_WRITE, and CG_MODE_MODIFY.
    ///
    /// returns CGNS file index number.
    pub fn open(path: &str, mode: CgMode) -> Result<c_int, Box<dyn Error>> {
        let mut i_file: c_int = 0;

        unsafe {
            check(cg_open(
                CString::new(path)?.as_ptr(),
                mode.value() as c_int,
                &mut i_file as *mut c_int,
            ))?;
        }

        Ok(i_file)
    }

    /// Close a CGNS file.
    ///
    /// Arguments:
    /// * `i_file` : CGNS file index number.
    pub fn close(i_file: &mut c_int) -> Result<(), Box<dyn Error>> {
        unsafe {
            check(cg_close(*i_file))?;
        }
        *i_file = 0;
        Ok(())
    }

    /// Create and/or write to a CGNS base node.
    ///
    /// Arguments:
    /// * `i_file` : CGNS file index number.
    /// * `base_name` : Name of the base.
    /// * `cell_dim` : Dimension of the cells; 3 for volume cells, 2 for surface cells and 1 for line cells.
    /// * `phys_dim` : Number of coordinates required to define a vector in the field.
    ///
    /// returns Base index number, where 1 ≤ B ≤ nbases.
    pub fn base_write(
        i_file: c_int,
        base_name: &str,
        cell_dim: c_int,
        phys_dim: c_int,
    ) -> Result<c_int, Box<dyn Error>> {
        check_str_len(base_name)?;

        let mut i_base: c_int = 0;
        unsafe {
            check(cg_base_write(
                i_file,
                CString::new(base_name)?.as_ptr(),
                cell_dim,
                phys_dim,
                &mut i_base,
            ))?;
        }

        Ok(i_base)
    }

    #[allow(dead_code)]
    #[derive(PartialEq)]
    pub enum ZoneType {
        Null,
        UserDefined,
        Structured,
        Unstructured,
    }

    impl ZoneType {
        fn value(&self) -> ZoneType_t {
            match self {
                ZoneType::Null => ZoneType_t_ZoneTypeNull,
                ZoneType::UserDefined => ZoneType_t_ZoneTypeUserDefined,
                ZoneType::Structured => ZoneType_t_Structured,
                ZoneType::Unstructured => ZoneType_t_Unstructured,
            }
        }
    }

    // TODO add CG_BUILD_64BIT detection
    // pub type CgSizeT = cgsize_t;
    pub type CgSizeT = c_int;

    /// Create and/or write to a zone node.
    ///
    /// Arguments:
    /// * `i_file` : CGNS file index number.
    /// * `i_base` : Base index number, where 1 ≤ B ≤ nbases.
    /// * `zone_name` : Name of the zone.
    /// * `size` : 	Number of vertices, cells, and boundary vertices in each (index)-dimension. For structured grids, the dimensions have unit stride in the array (e.g., [NVertexI, NVertexJ, NVertexK, NCellI, NCellJ, NCellK, NBoundVertexI, NBoundVertexJ, NBoundVertexK]).
    ///             Note that for unstructured grids, the number of cells is the number of highest order elements. Thus, in three dimensions it's the number of 3-D cells, and in two dimensions it's the number of 2-D cells.
    ///             Also for unstructured grids, if the nodes are sorted between internal nodes and boundary nodes, the optional parameter NBoundVertex must be set equal to the number of boundary nodes. By default, NBoundVertex equals zero, meaning that the nodes are unsorted.
    ///             Note that a non-zero value for NBoundVertex only applies to unstructured grids. For structured grids, the NBoundVertex parameter always equals 0 in all directions.
    ///             Mesh Type         Size
    ///             3D structured     NVertexI, NVertexJ, NVertexK
    ///                               NCellI, NCellJ, NCellK
    ///                               NBoundVertexI = 0, NBoundVertexJ = 0, NBoundVertexK = 0
    ///             2D structured     NVertexI, NVertexJ
    ///                               NCellI, NCellJ
    ///                               NBoundVertexI = 0, NBoundVertexJ = 0
    ///             3D unstructured   NVertex, NCell3D, NBoundVertex
    ///             2D unstructured	  NVertex, NCell2D, NBoundVertex
    /// * `type_` : 	Type of the zone. The admissible types are Structured and Unstructured.
    ///
    /// returns Zone index number, where 1 ≤ Z ≤ nzones.
    pub fn zone_write(
        i_file: c_int,
        i_base: c_int,
        zone_name: &str,
        size: &Vec<CgSizeT>,
        type_: &ZoneType,
    ) -> Result<c_int, Box<dyn Error>> {
        assert!(*type_ == ZoneType::Structured || *type_ == ZoneType::Unstructured);

        let mut i_zone: c_int = 0;

        unsafe {
            check(cg_zone_write(
                i_file,
                i_base,
                CString::new(zone_name)?.as_ptr(),
                size.as_ptr(),
                type_.value(),
                &mut i_zone,
            ))?;
        }

        Ok(i_zone)
    }

    #[allow(dead_code)]
    #[derive(PartialEq)]
    enum DataType {
        Null,
        UserDefined,
        Integer,
        RealSingle,
        RealDouble,
        Character,
        LongInteger,
    }

    impl DataType {
        fn value(&self) -> u32 {
            match self {
                DataType::Null => DataType_t_DataTypeNull,
                DataType::UserDefined => DataType_t_DataTypeUserDefined,
                DataType::Integer => DataType_t_Integer,
                DataType::RealSingle => DataType_t_RealSingle,
                DataType::RealDouble => DataType_t_RealDouble,
                DataType::Character => DataType_t_Character,
                DataType::LongInteger => DataType_t_LongInteger,
            }
        }
    }

    /// Write grid coordinates.
    ///
    /// Arguments:
    /// * `i_file` : CGNS file index number.
    /// * `i_base` : Base index number, where 1 ≤ B ≤ nbases.
    /// * `i_zone` : Zone index number, where 1 ≤ Z ≤ nzones.
    /// * `data_type` : Data type of the coordinate array written to the file. Admissible data types for a coordinate array are RealSingle and RealDouble.
    /// * `coord_name` : Name of the coordinate array. It is strongly advised to use the SIDS nomenclature conventions when naming the coordinate arrays to insure file compatibility.
    /// * `coord_array` : Array of coordinate values.
    ///
    /// returns Coordinate array index number, where 1 ≤ C ≤ ncoords.
    pub fn coord_write(
        i_file: c_int,
        i_base: c_int,
        i_zone: c_int,
        coord_name: &str,
        coord_array: &[f64],
    ) -> Result<c_int, Box<dyn Error>> {
        // assert!(*data_type == DataType::RealSingle || *data_type ==
        // DataType::RealDouble);

        // TODO extend to allow c_float

        let data_type = DataType::RealDouble;

        let mut i_coord: c_int = 0;
        unsafe {
            check(cg_coord_write(
                i_file,
                i_base,
                i_zone,
                data_type.value(),
                CString::new(coord_name)?.as_ptr(),
                coord_array.as_ptr() as *const c_void,
                &mut i_coord,
            ))?;
        }

        Ok(i_coord)
    }
}

// /// cgns file
// pub struct File {
//     /// file handle
//     i_file: c_int,
// }

// impl File {
//     pub fn new(path: &str) -> Result<Self, Box<dyn Error>> {
//         Self::_open(path, CgMode::Write)
//     }

//     /// open file to modify
//     fn open(path: &str) -> Result<Self, Box<dyn Error>> {
//         Self::_open(path, CgMode::Modify)
//     }

//     /// wrapper for cg_close. Called in drop method.
//     fn _close(&self) {
//         unsafe {
//             let ier = cg_close(self.i_file);
//         }
//     }
// }

// impl Drop for File {
//     fn drop(&mut self) {
//         self._close();
//     }
// }
