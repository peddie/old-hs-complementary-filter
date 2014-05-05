{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}

{-|
Module      :  AttCmpl
Copyright   :  (c) Matthew Peddie 2013
License     :  Public Domain

Maintainer  :  peddie@alum.mit.edu
Stability   :  experimental
Portability :  GHC

Attitude complementary filter.

-}

module AttCmpl where

import Linear
import Data.List (foldl')

import Data.Random

-- | The filter state consists of a quaternion @q_r2b@ describing the
-- rotation from the reference frame (e.g. NED, LVLH) to the body
-- frame and a 3-vector @β@ of gyro biases.
type FilterState a = (Quaternion a, V3 a)

type SystemState a = (Quaternion a, V3 a)

-- | Normalize a quaternion.
normalize' :: Floating a => Quaternion a -> Quaternion a
normalize' q = fmap (* normInv) q
  where
    normInv = 1 / norm q

-- | Integrate a quaternion @q@ given a rotation vector @ω@.
qint :: (Floating a, Fractional a) => a  -- ^ Δt
     -> Quaternion a  -- ^ q
     -> V3 a          -- ^ ω
     -> Quaternion a
qint dt (Quaternion q0 (V3 q1 q2 q3)) (V3 x y z) = normalize' $ Quaternion q0' (V3 q1' q2' q3')
  where
    q0' = q0 + dt*0.5*(-x*q1 - y*q2 - z*q3)
    q1' = q1 + dt*0.5*( x*q0 + z*q2 - y*q3)
    q2' = q2 + dt*0.5*( y*q0 - z*q1 + x*q3)
    q3' = q3 + dt*0.5*( z*q0 + y*q1 - x*q2)

-- | Quaternion complementary filter step.  Given a set of measured
-- vectors and reference vectors with gains, a measured angular
-- velocity vector, some filter gains and a filter state, it outputs
-- the new filter state.
qcomplStep :: (Conjugate a, RealFloat a) =>
              a                     -- ^ Δt
           -> [(V3 a, V3 a, a)]  -- ^ List of triples: (measured vector (body frame), reference vector (reference frame), gain)
           -> V3 a               -- ^ Gyro measurement
           -> (a, a)             -- ^ Filter gains (kp, ki)
           -> FilterState a      -- ^ Initial state (q, β)
           -> FilterState a
qcomplStep dt lowfreq gyro (kp, ki) (oldq, oldbeta) = (q, beta)
  where
    (gainsum, lf) = foldl' vecerr (0, V3 0 0 0) lowfreq
    vecerr (kacc, vacc) (vec, refvec, gain) = (kacc + gain, vacc - (gain *^ vec `cross` rotate oldq refvec))
    lferr = lf ^* (kp / gainsum)
    eq = lferr + gyro - oldbeta
    q = qint dt oldq eq
    bup = -ki *^ lferr
    beta = oldbeta + bup ^* dt

normalErrQuat :: (Floating a, Distribution Normal a) => a -> RVarT m (Quaternion a)
normalErrQuat sig = normalErrVec sig >>= return . normalize' . Quaternion (1 - sig**2)

normalErrVec :: (Floating a, Distribution Normal a) => a -> RVarT m (V3 a)
normalErrVec sig = do
  x <- stdNormalT
  y <- stdNormalT
  z <- stdNormalT
  let evec = V3 x y z
  return $ evec ^* (sig / norm evec)

runFilterConst :: (Enum a, RealFloat a, Conjugate a, Distribution Normal a) =>
                  [(V3 a, a, a)]
               -> (V3 a, a, a)
               -> a
               -> a
               -> SystemState a
               -> a
               -> RVar [(a, FilterState a)]
runFilterConst lfvecs (omega, ki, sigg) kp dt (qtrue, betatrue) tmax = do
  meas <- mapM (\(v, k, sig) -> normalErrQuat sig >>= \qerr ->
                 return (rotate (qtrue * qerr) v, v, k)) lfvecs
  gyronoise <- normalErrVec sigg
  let upd = qcomplStep dt meas (omega + gyronoise + betatrue) (kp, ki)
      ts = [0,dt..tmax]
      traj = foldl' (\st@((_, prev):_) t -> (t, upd prev) : st) [(0, (Quaternion 1 (V3 0 0 0), V3 0 0 0))] ts
  return traj
      -- meas = map (\(v, k) -> (rotate qtrue v, v, k)) lfvecs
