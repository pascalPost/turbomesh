const std = @import("std");

const GeomTag = enum { spline, point, line };

const Geom = union(GeomTag) {};

pub const MeshGenerator = struct {
    // /// tables (from csv), single values, ...
    // data: [] ,
    //
    // /// Points, Lines,
    // geometry: [],
    //
    // /// edges, ...
    // discrete: []

    // fn addSpline() void
};
