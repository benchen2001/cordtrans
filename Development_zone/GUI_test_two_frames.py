# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 11:29:43 2015

@author: user
"""

from Tkinter import *

root = Tk()
#%% Frame 1
Frame1 = Frame(bg="blue",width=176, height=157)
Frame1.grid(row=0,column=0)
la1f1=Label(Frame1,text="test1")
la1f1.grid(row=0,column=0)
la2f1=Label(Frame1,text="test2")
la2f1.grid(row=1,column=1)


#%% Frame 2
Frame2 = Frame(root,bg="green",width=276, height=257)
Frame2.grid(row=0,column=1)
bt_TD1 = Label(Frame2, text="Demo - precesion & nutation - regular motion")#, command=TD1callback)
bt_TD1.pack()

root.mainloop()