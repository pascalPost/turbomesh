use crate::types::Scalar;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec2d(pub Scalar, pub Scalar);

impl Vec2d {
    pub fn abs(&self) -> Scalar {
        (self.0 * self.0 + self.1 * self.1).sqrt()
    }
}

impl std::ops::Add for Vec2d {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0, self.1 + rhs.1)
    }
}

impl std::ops::Sub for Vec2d {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0, self.1 - rhs.1)
    }
}

impl std::ops::Mul<Vec2d> for Scalar {
    type Output = Vec2d;

    fn mul(self, rhs: Vec2d) -> Self::Output {
        Vec2d(self * rhs.0, self * rhs.1)
    }
}
