#!/bin/python3

import turtle
from tkinter import *

wn = turtle.Screen()
wn.setup(400,600)
wn.bgcolor("white")
s = turtle.Turtle()
s.hideturtle()
s.speed(0)

for i in range(4):
    s.begin_fill()
    s.forward(100)
    s.left(90)
    s.pencolor("pink")
    s.forward(100)
    s.fillcolor("pink")
    s.end_fill()

for i in range(4):
    s.begin_fill()
    s.forward(90)
    s.left(90)
    s.pencolor("purple")
    s.forward(90)
    s.fillcolor("purple")
    s.end_fill()

for i in range(4):
    s.begin_fill()
    s.forward(80)
    s.left(90)
    s.pencolor("green")
    s.forward(80)
    s.fillcolor("green")
    s.end_fill()

for i in range(4):
    s.begin_fill()
    s.forward(70)
    s.left(90)
    s.pencolor('orange')
    s.forward(70)
    s.fillcolor("orange")
    s.end_fill()

for i in range(4):
    s.begin_fill()
    s.forward(60)
    s.left(90)
    s.pencolor("yellow")
    s.forward(60)
    s.fillcolor("yellow")
    s.end_fill()
        
for i in range(4):
    s.begin_fill()
    s.forward(50)
    s.left(90)
    s.pencolor("grey")
    s.forward(50)
    s.fillcolor("grey")
    s.end_fill()

for i in range(4):
    s.penup()
    s.goto(0, 180)
    s.write("ensembleamppred", move=True, align="center", font=("Courier", 15, "normal"))
    s.goto(0, 160)
    s.write("macrel", move=True, align="center", font=("Courier", 15, "normal"))
    s.goto(0, 140)
    s.write("amplify", move=True, align="center", font=("Courier", 15, "normal"))
    s.goto(0, 120)
    s.write("hmmsearch", move=True, align="center", font=("Courier", 15, "normal"))
    s.goto(0, 100)
    s.write("ampir", move=True, align="center", font=("Courier", 15, "normal"))
    s.goto(0, 80)
    s.write("neubi", move=True, align="center", font=("Courier", 15, "normal"))
    s.color("Black")

for i in range(4):
    s.penup()
    s.goto(0, -40)
    s.write("AMP_COMBI",font=("Arial", 25, "bold"),move=True,align="center")
    s.pencolor('Black')
    s.color('Black')

s.getscreen().getcanvas().postscript(file='outputname.ps')

#for i in range(4):
    #s.forward(60)
    #s.left(60)
    #s.penup()
    #s.home()
    #s.dot()
    #s.write("AMP_combi",font=("Helvetica", 25, "normal"),align="center")
    #s.color("Black")
#turtle.write("AMP_combi",font=("Helvetica", 25, "normal"),align="center")
#turtle.color('deep pink')
#style = ('Courier', 30, 'italic')
#turtle.write('Hello!', font=style, align='center')
#turtle.hideturtle()
turtle.done()