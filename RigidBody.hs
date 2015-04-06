{-# OPTIONS_GHC -Wall #-}

module RigidBody where

import Numeric.LinearAlgebra.HMatrix as H
import Numeric.GSL.ODE

{- Quaternions, embedded in HMatrix vectors -}

qId = quat 1 0 0 0

quatInv v = fromList [q0, -q1, -q2, -q3]
  where
    q0 = v ! 0
    q1 = v ! 1
    q2 = v ! 2
    q3 = v ! 3
{-# INLINE quatInv #-}

quatNorm = norm_2
{-# INLINE quatNorm #-}

quatNormalize v = (1.0 / quatNorm v) `scale` v
{-# INLINE quatNormalize #-}

quat q0 q1 q2 q3 = quatNormalize $ fromList [q0, q1, q2, q3]

-- | Quaternion multiplication.
mul qa qb
  | q0 < 0 = quat (-q0) (-q1) (-q2) (-q3)
  | otherwise = quat q0 q1 q2 q3
  where
    qa0 = qa ! 0
    qa1 = qa ! 1
    qa2 = qa ! 2
    qa3 = qa ! 3
    qb0 = qb ! 0
    qb1 = qb ! 1
    qb2 = qb ! 2
    qb3 = qb ! 3
    q0 = qa0*qb0 - qa1*qb1 - qa2*qb2 - qa3*qb3;
    q1 = qa0*qb1 + qa1*qb0 + qa2*qb3 - qa3*qb2;
    q2 = qa0*qb2 - qa1*qb3 + qa2*qb0 + qa3*qb1;
    q3 = qa0*qb3 + qa1*qb2 - qa2*qb1 + qa3*qb0;

{- Quaternion kinematics -}

-- | Create a linearized "differential quaternion" (for use with an
-- external integrator)
dq q w = fromList [q0', q1', q2', q3']
  where
    q0' = 0.5*(-wx*q1 - wy*q2 - wz*q3)
    q1' = 0.5*( wx*q0 + wz*q2 - wy*q3)
    q2' = 0.5*( wy*q0 - wz*q1 + wx*q3)
    q3' = 0.5*( wz*q0 + wy*q1 - wx*q2)
    q0 = q ! 0
    q1 = q ! 1
    q2 = q ! 2
    q3 = q ! 3
    wx = w ! 0
    wy = w ! 1
    wz = w ! 2

-- | Do a direct Euler step.
qEuler dt q w = quat q0' q1' q2' q3'
  where
    q0' = q0*c    - q1*wx*s - q2*wy*s - q3*wz*s;
    q1' = q0*wx*s + q1*c    + q2*wz*s - q3*wy*s;
    q2' = q0*wy*s - q1*wz*s + q2*c    + q3*wx*s;
    q3' = q0*wz*s + q1*wy*s - q2*wx*s + q3*c;

    wnorm = sqrt $ w `dot` w
    wdt2 = 0.5 * wnorm * dt
    c = cos wdt2
    s = sin wdt2 / wnorm

    q0 = q ! 0
    q1 = q ! 1
    q2 = q ! 2
    q3 = q ! 3
    wx = w ! 0
    wy = w ! 1
    wz = w ! 2

{- Spatial rotation with quaternions -}

quatToDCM q = (3><3) [d0, d1, d2, d3, d4, d5, d6, d7, d8]
  where
    d0 = q0*q0 + q1*q1 - q2*q2 - q3*q3
    d3 = 2*(q1*q2 - q0*q3)
    d6 = 2*(q1*q3 + q0*q2)

    d1 = 2*(q1*q2 + q0*q3)
    d4 = q0*q0 - q1*q1 + q2*q2 - q3*q3
    d7 = 2*(q2*q3 - q0*q1)

    d2 = 2*(q1*q3 - q0*q2)
    d5 = 2*(q2*q3 + q0*q1)
    d8 = q0*q0 - q1*q1 - q2*q2 + q3*q3
    q0 = q ! 0
    q1 = q ! 1
    q2 = q ! 2
    q3 = q ! 3

rotVecByQuat v q = quatToDCM q #> v

rotVecByQuatInv v q = tr (quatToDCM q) #> v

{- Rigid-body dynamics

'j' is moment of inertia, typically written I or J.

'q_i2b' means "the quaternion that represents the rotation from the
inertial reference frame to the body frame."

'w_bi_b' means "the angular velocity vector describing the rotation of
the body frame with respect to the inertial frame, expressed in the
body frame."

-}

rigidBodyStep (j_b, j_b_inv) (torque_i, torque_b) _ state =
  vjoin [q'_i2b, a_bi_b, power, impulse]
  where
    [q_i2b, w_bi_b] = takesV [4, 3] state
    -- rotate inertial-frame torques into the body frame
    torque = torque_b + (torque_i `rotVecByQuat` q_i2b)
    -- Compute rotational acceleration in the body frame
    a_bi_b = j_b_inv #> (torque - (w_bi_b `cross` (j_b #> w_bi_b)))
    -- Compute differential quaternion
    q'_i2b = dq q_i2b w_bi_b
    -- Compute power; this lets us track total energy, which we can
    -- use as a check that our dynamics are behaving sensibly.
    power = scalar $ w_bi_b `dot` torque
    -- Compute impulse; this lets us track angular momentum, which
    -- serves as another check.
    impulse = scalar . norm_2 $ j_b #> a_bi_b
    -- TODO(MP): Test this under a wider range of conditions.
    --
    -- Compute a divergence term to correct GSL state updates.  (This
    -- has no physical meaning.)
    -- q_i2b = quatNormalize q_i2b_pre
    -- qdiv = 1e-6 * (q_i2b - q_i2b_pre)

odeV :: (Double -> Vector Double -> Vector Double)
     -> Vector Double -> Vector Double -> Matrix Double
odeV = odeSolveV MSAdams 0.001 1e-9 1e-9

testJ = (j_b, inv j_b)
  where
    j_b = ident 3

testTorque :: (Vector Double, Vector Double)
testTorque = (konst 0 3, konst 0 3)

testQ = fromList [1, 0, 0, 0]
testW = fromList [1, 0, 0]
testEnergy = scalar $ 0.5 * (norm_2 $ fst testJ #> testW)**2
testMomentumNorm = scalar . norm_2 $ fst testJ #> testW

testState = vjoin [testQ, testW, testEnergy, testMomentumNorm]
