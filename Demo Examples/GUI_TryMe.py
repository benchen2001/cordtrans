# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 16:02:15 2015

@author: whymranderson.blogspot.tw
"""

import Tkinter as TK
#from PIL import ImageTk, Image

master = TK.Tk()
master.wm_title("Gyroscope Motion, Rigid Body Intergrator, and Orientation Estimation Platform")
master.minsize(width=555, height=150)


def TD1callback():
    variables= {'testparameter':e1.get()}
    execfile( "Gyroscope-TeachDemo-1_Gyro_Cube_Motion_A.py", variables)#    variables2 
    print variables('line_ani')
def TD2callback():
    variables= {}
    execfile( "Gyroscope_SpaceBodyCone.py", variables)#    variables2 )
    print variables('line_ani')
def TD3callback():
    variables= {}
    execfile( "Gyroscope-TeachDemo-1_Gyro_Cube_Motion_circular.py", variables)#    variables2 )
    print variables('line_ani')
def TD4callback():
    variables= {}
    execfile( "Gyroscope-TeachDemo-3.py", variables)#    variables2 )
    print variables('line_ani')    
def TD5callback():
    variables= {}
    execfile( "Gyroscope-TeachDemo-4.py", variables)#    variables2 )
    print variables('line_ani')
def TD6callback():
    variables= {}
    execfile( "Gyroscope-TeachDemo-5-NoiseIncludeInCMethod.py", variables)#    variables2 )
    print variables('line_ani')    
def TD7callback():
    variables= {}
    execfile( "Gyroscope-TeachDemo-1_Gyro_Cube_Motion_ring.py", variables)#    variables2 )
    print variables('line_ani')   
def TD8callback():
    variables= {}
    execfile( "Gyroscope-TeachDemo-1_Gyro_Cube_Motion_wave.py", variables)#    variables2 )
    print variables('line_ani')
def TD9callback():
    variables= {}
    execfile( "Gyroscope-TeachDemo-5-NoiseIncludeInCMethod_still.py", variables)#    variables2 )
    print variables('line_ani')
def TD10callback():
    variables= {}
    execfile( "AngularVelocityTrail_in_body_frame_compare_ABmethods.py", variables)#    variables2 )
    print variables('line_ani')


bt_TD1 = TK.Button(master, text="Demo - precesion & nutation - regular motion", command=TD1callback)
bt_TD1.pack()
bt_TD7 = TK.Button(master, text="Demo - precesion & nutation - classical ring motion", command=TD7callback)
bt_TD7.pack()
bt_TD8 = TK.Button(master, text="Demo - precesion & nutation - classical wave motion", command=TD8callback)
bt_TD8.pack()
bt_TD2 = TK.Button(master, text="Demo(static) - Space Body Cone", command=TD2callback)
bt_TD2.pack()
bt_TD3 = TK.Button(master, text="Demo - circular locus", command=TD3callback)
bt_TD3.pack()
bt_TD4 = TK.Button(master, text="Demo - AB methods compare", command=TD4callback)
bt_TD4.pack()
bt_TD5 = TK.Button(master, text="Demo - BC methods compare", command=TD5callback)
bt_TD5.pack()
bt_TD6 = TK.Button(master, text="Demo - Noise included (moving)", command=TD6callback)
bt_TD6.pack()
bt_TD9 = TK.Button(master, text="Demo - Noise included (still)", command=TD9callback)
bt_TD9.pack()
bt_TD10 = TK.Button(master, text="Demo - Angular Velocity Trail in the body frame", command=TD10callback)
bt_TD10.pack()

#%% Set initial conditions
labelte1=TK.StringVar()
labelte1.set("enter wx here")
labele1=TK.Label(master, textvariable=labelte1)
labele1.pack()
e1 = TK.Entry(master)
e1.pack()

#%% Text
text = TK.Text(master,height=2,font=("Helvetica", 10, "italic"))
text.insert(TK.INSERT, "\nProduced by whymranderson.blogspot.tw")
#text.insert(TK.END, "Bye Bye.....")
text.pack()

#%% Image
'''
path = 'C:/Documents and Settings/user/My Documents/tony/Scripts/GyroDocs/new_gyro_fig.png'

img = ImageTk.PhotoImage(Image.open(path))
panel = TK.Label(master, image = img)
panel.pack(side = "bottom", fill = "both", expand = "yes")
'''

TK.mainloop()


