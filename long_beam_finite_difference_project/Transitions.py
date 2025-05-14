# -*- coding: utf-8 -*-
"""
Created on Sun May 21 18:22:48 2023

@author: mbluu
"""

# -*- coding: utf-8 -*-
"""
Created on Mon May 15 17:40:36 2023

@author: mbluu
"""

from manim import *

INITIAL_SHIFT = 2


def Develop_Img(name,DIR):
    img = ImageMobject(name)
    img.scale(0.4)
    img.to_corner(DIR)
    
    return img

class Transition_MovingCamera(MovingCameraScene):
    def construct(self):
        img1 = Develop_Img('Beam_Derivation_KeyFrame.png',UP+LEFT)
        box1 = SurroundingRectangle(img1, corner_radius=0.2)
        group1 = Group(img1,box1)
        
        img11 = Develop_Img('Beam_Derivation_KeyFrame2.png',UP+LEFT)
        box11 = SurroundingRectangle(img11, corner_radius=0.2)
        group11 = Group(img11,box11)
        
        img2 = Develop_Img('FE_Difference_KeyFrame.png',UP+RIGHT)
        box2 = SurroundingRectangle(img2, corner_radius=0.2)
        group2 = Group(img2,box2)

        img3 = Develop_Img('FE_Difference_Inverse_KeyFrame.png',DOWN+LEFT)
        box3 = SurroundingRectangle(img3, corner_radius=0.2)
        group3 = Group(img3,box3)

        img4 = Develop_Img('Wave_Prop_KeyFrame.png',DOWN+RIGHT)
        box4 = SurroundingRectangle(img4, corner_radius=0.2)
        group4 = Group(img4,box4)

        pause = 3
        self.play(FadeIn(img1),FadeIn(img2),FadeIn(img3),FadeIn(img4))
        self.wait(1)
        self.play(Create(box1),Create(box2),Create(box3),Create(box4))
        self.wait(pause)
        self.wait(1)

        self.camera.frame.save_state()
        self.play(self.camera.frame.animate.set(width=group1.width * 1.1).move_to(group1))
        self.wait(pause)
        self.play(Restore(self.camera.frame))
        self.wait(1)

        self.play(self.camera.frame.animate.set(width=group3.width * 1.1).move_to(group2))
        self.wait(pause)    
        self.play(Restore(self.camera.frame))
        self.wait(1)
        
        self.play(self.camera.frame.animate.set(width=group3.width * 1.1).move_to(group3))
        self.wait(pause)    
        self.play(Restore(self.camera.frame))
        self.wait(1)

        self.play(self.camera.frame.animate.set(width=group2.width * 1.1).move_to(group4))
        self.wait(pause)    
        self.play(Restore(self.camera.frame),FadeOut(group1),FadeIn(group11))
        self.wait(pause)    

        #self.play(FadeOut(group1),FadeIn(group11))
        self.play(self.camera.frame.animate.set(width=group11.width * 1.0).move_to(group11,UP),FadeOut(box11),FadeOut(box1))
