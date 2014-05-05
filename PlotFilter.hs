{-# OPTIONS_GHC -Wall #-}

{-|
Module      :  PlotFilter
Copyright   :  (c) Matthew Peddie 2013
License     :  Public Domain

Maintainer  :  peddie@alum.mit.edu
Stability   :  experimental
Portability :  GHC

Test driver for the complementary filter

-}

module Main where

import Graphics.Rendering.Chart
import Data.Colour
import Data.Colour.Names
import Data.Default.Class
import Graphics.Rendering.Chart.Gtk
import Control.Lens

import AttCmpl
import Linear

import Data.Random
import Data.Random.Source.PureMT
import Data.IORef

sim :: RVar [(Double, FilterState Double)]
sim = runFilterConst
       [(V3 1 0 0, 1, 0.05), (V3 0 1 0, 1, 0.1)]
       (V3 0 0 0, 0.01, 0.18) 5 0.2
       (Quaternion 0.8 (V3 0 0.2 0), V3 0.22 (-0.05) 0.1)
       500

setLinesBlue :: PlotLines a b -> PlotLines a b
setLinesBlue = plot_lines_style  . line_color .~ opaque blue

biasChart :: [(Double, FilterState Double)] -> Renderable ()
biasChart traj = toRenderable layout
  where
    biasx = plot_lines_values .~ [traj & mapped._2 %~ (^._2._x)]
            $ plot_lines_style  . line_color .~ opaque blue
            $ plot_lines_title .~ "bias x"
            $ def
    biasy = plot_lines_values .~ [traj & mapped._2 %~ (^._2._y)]
            $ plot_lines_style  . line_color .~ opaque green
            $ plot_lines_title .~ "bias y"
            $ def
    biasz = plot_lines_values .~ [traj & mapped._2 %~ (^._2._z)]
            $ plot_lines_style  . line_color .~ opaque red
            $ plot_lines_title .~ "bias z"
            $ def

    layout = layout_title .~ "Bias estimates vs. time"
           $ layout_plots .~ [toPlot biasx, toPlot biasy, toPlot biasz]
           $ def

quatChart :: [(Double, FilterState Double)] -> Renderable ()
quatChart traj = toRenderable layout
  where
    q0 = plot_lines_values .~ [traj & mapped._2 %~ (^._1._e)]
         $ plot_lines_style  . line_color .~ opaque pink
         $ plot_lines_title .~ "q0"
         $ def
    q1 = plot_lines_values .~ [traj & mapped._2 %~ (^._1._i)]
         $ plot_lines_style  . line_color .~ opaque blue
         $ plot_lines_title .~ "q1"
         $ def
    q2 = plot_lines_values .~ [traj & mapped._2 %~ (^._1._j)]
         $ plot_lines_style  . line_color .~ opaque green
         $ plot_lines_title .~ "q2"
         $ def
    q3 = plot_lines_values .~ [traj & mapped._2 %~ (^._1._k)]
         $ plot_lines_style  . line_color .~ opaque red
         $ plot_lines_title .~ "q3"
         $ def

    layout = layout_title .~ "Quaternion parameters vs. time"
           $ layout_plots .~ map toPlot [q0, q1, q2, q3]
           $ def

main :: IO ()
main = do
  src <- newIORef =<< newPureMT
  traj <- runRVar sim src
  renderableToWindow (biasChart traj) 800 600
  renderableToWindow (quatChart traj) 800 600
