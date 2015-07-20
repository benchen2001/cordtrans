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

Frame1 = TK.Frame(master,bg="blue")#,width=176, height=157)
Frame1.grid(row=0,column=0)
Frame2 = TK.Frame(master,bg="green")#,width=176, height=157)
Frame2.grid(row=0,column=1)



#%% callbacks
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
def fetch_n_execute_callback():
    variables= {'set_wx':E1.get(),'set_wy':E2.get(),'set_g':E3.get(),'set_tf':E4.get()}
    execfile( "Gyroscope-TeachDemo-custumed_parameters.py", variables)#    variables2 )
    print variables('line_ani')

#%% Frame 2 widgets

f2title = TK.Label(Frame2,text="Demo Examples",bg="white")
f2title.pack()
bt_TD1 = TK.Button(Frame2, text="Demo - precesion & nutation - regular motion", command=TD1callback)
bt_TD1.pack()
bt_TD7 = TK.Button(Frame2, text="Demo - precesion & nutation - classical ring motion", command=TD7callback)
bt_TD7.pack()
bt_TD8 = TK.Button(Frame2, text="Demo - precesion & nutation - classical wave motion", command=TD8callback)
bt_TD8.pack()
bt_TD2 = TK.Button(Frame2, text="Demo(static) - Space Body Cone", command=TD2callback)
bt_TD2.pack()
bt_TD3 = TK.Button(Frame2, text="Demo - circular locus", command=TD3callback)
bt_TD3.pack()
bt_TD4 = TK.Button(Frame2, text="Demo - AB methods compare", command=TD4callback)
bt_TD4.pack()
bt_TD5 = TK.Button(Frame2, text="Demo - BC methods compare", command=TD5callback)
bt_TD5.pack()
bt_TD6 = TK.Button(Frame2, text="Demo - Noise included (moving)", command=TD6callback)
bt_TD6.pack()
bt_TD9 = TK.Button(Frame2, text="Demo - Noise included (still)", command=TD9callback)
bt_TD9.pack()
bt_TD10 = TK.Button(Frame2, text="Demo - Angular Velocity Trail in the body frame", command=TD10callback)
bt_TD10.pack()

#%% Frame 1 Set initial conditions
f1title = TK.Label(Frame1,text="set parameters",bg="white")
f1title.grid(row=0,column=0,columnspan=2)
TK.Label(Frame1,text="wx = (radian/second)").grid(row=1,column=0)
E1=TK.Entry(Frame1)
E1.grid(row=1,column=1)
TK.Label(Frame1,text="wy=").grid(row=2,column=0)
E2 = TK.Entry(Frame1)
E2.grid(row=2,column=1)
TK.Label(Frame1,text="gravity").grid(row=3,column=0)
E3 = TK.Entry(Frame1)
E3.grid(row=3,column=1)
TK.Label(Frame1,text="stop time (in seconds) =").grid(row=4,column=0)
E4 = TK.Entry(Frame1)
E4.grid(row=4,column=1)

B2=TK.Button(Frame1,text="run",command=fetch_n_execute_callback)
B2.grid(row=5,columnspan=2)
#E2.get()

#%% Text
text = TK.Text(master,height=2,font=("Helvetica", 10, "italic"))
text.insert(TK.INSERT, "\nProduced by whymranderson.blogspot.tw")
#text.insert(TK.END, "Bye Bye.....")
text.grid(row=1,column=0,columnspan=2)

#%% Image
'''
path = 'C:/Documents and Settings/user/My Documents/tony/Scripts/GyroDocs/new_gyro_fig.png'

img = ImageTk.PhotoImage(Image.open(path))
panel = TK.Label(master, image = img)
panel.pack(side = "bottom", fill = "both", expand = "yes")
'''

TK.mainloop()


