# coding: utf-8

q=linspace(0,1,0.1)
get_ipython().magic('pylab')
s,t=ix_(t,t.T)
s,t=ix_(q,q.T)
q=linspace(0,1,0.1)
s,t=ix_(q,q.T)
s
t
s,t=meshgrid(q,q)
s
q
q=linspace(0,1,11)
q
s,t=meshgrid(q,q)
s
figure(1)
plot(x,y)
plot(s,t)
plot(s.T,t.T)
plot(s.T,t.T,'ok')
clf()
plot(s.T,t.T,'ok')
plot(s.T,t.T,'k-')
cls()
clf()
plot(s.T,t.T,'k-')
plot(s,t,'k-')
figure(2)
def ft(x,y,a,b,c):
    x = (s+1)/4 * ((1-t)*a + (t+1)*(a+c))
    y = (t+1)/4 * ((s-1)*a + (s+1)*(c-b))
    return x,y
def ft(s,t,a,b,c):
    x = (s+1)/4 * ((1-t)*a + (t+1)*(a+c))
    y = (t+1)/4 * ((s-1)*a + (s+1)*(c-b))
    return x,y
x,y=ft(s,t,2,1,5)
plot(x,y,'k-')
cls()
clf()
plot(s,t,'k-')
plot(s.T,t.T,'k-')
def g(s,t):
    plot(s,t,'-k')
    plot(s.T,t.T,'-ok')
    
clf()
g(s,t)
figure(2)
g(x,y)
q=linspace(-1,1,21)
s,t=meshgrid(q,q)
figure(!)
clf()
figure(1)
clf()
g(s,t)
clf()
q=linspace(-1,1,11)
s,t=meshgrid(q,q)
g(s,t)
figure(2)
x,y=ft(s,t,2,1,5)
g(x,y)
def g(s,t):
    plot(s,t,'-k')
    plot(s.T,t.T,'-k')
    scatter(s,t,c=s+t)
    
figure(1)
clf()
g(s,t)
w=linspace(0,1,11)
u,v=meshgrid(w,w)
figure(3)
g(u,v)
figure(2)
clf()
def ft_uv(u,v,a,b,c)
def ft_uv(u,v,a,b,c):
    x = u*((1-v)*a+v*(a+c))
    y = v*((u-1)*a+u*(c-b))
    return x,y
g(ft_uv(u,v,2,1,5))
g(*ft_uv(u,v,2,1,5))
ft_uv(0,0,2,1,5)
ft_uv(1,0,2,1,5)
ft_uv(0,1,2,1,5)
def ft_uv(u,v,a,b,c):
    x = a*u*(1-v) + a*v*(u-1)
    y = b*v*(1-u) + u*(b+v*(c-b))
    return x,y
clf()
g(*ft_uv(u,v,2,1,5))
clf();g(*ft_uv(u,v,2,1,5))
clf();g(*ft_uv(u,v,1,1,5))
def g(s,t):
    plot(s,t,'-k')
    plot(s.T,t.T,'-k')
    
clf();g(*ft_uv(u,v,1,1,5))
clf();g(*ft_uv(u,v,1,3,5))
clf();g(*ft_uv(u,v,1,4,5))
clf();g(*ft_uv(u,v,1,0.5,6))
clf();g(*ft_uv(u,v,12,0.5,6))
clf();g(*ft_uv(u,v,12,3,6))
def ft(s,t,a,b,c):
    x = a/4*((s+1)*(1-t)+(t+1)*(1-s))
    y = b/4*((t+1)*(1-s)+(s+1)(2+(t+1)*(c/b-1))
    return x,y)
def ft(s,t,a,b,c):
    x = a/4*((s+1)*(1-t)+(t+1)*(1-s))
    y = b/4*((t+1)*(1-s)+(s+1)(2+(t+1)*(c/b-1)))
    return x,y
figure(4);g(*ft(s,t,2,1,5))
def ft(s,t,a,b,c):
    x = a/4*((s+1)*(1-t)+(t+1)*(1-s))
    y = b/4*((t+1)*(1-s)+(s+1)*(2+(t+1)*(c/b-1)))
    return x,y
clf();g(*ft(s,t,12,3,6))
clf();g(*ft(s,t,12,3,6))
s
t
def ft(s,t,a,b,c):
    x = a/4*((s+1)*(1-t)+(t+1)*(s-1))
    y = b/4*((t+1)*(1-s)+(s+1)*(2+(t+1)*(c/b-1)))
    return x,y
clf();g(*ft(s,t,12,3,6))
clf();g(*ft(s,t,12,1,6))
clf();g(*ft(s,t,10,1,6))
clf();g(*ft(s,t,10,1,5))
import sympy as sp
ss,tt=sp.var('s,t')
ss
tt
a,b,c=sp.var('a,b,c')
sp.diff(ft(ss,tt,a,b,c)[0],ss)
sp.diff(ft(ss,tt,a,b,c)[0],tt)
sp.diff(ft(ss,tt,a,b,c)[1],ss)
sp.diff(ft(ss,tt,a,b,c)[1],tt)
sp.simplify(sp.diff(ft(ss,tt,a,b,c)[1],ss))
sp.simplify(sp.diff(ft(ss,tt,a,b,c)[1],tt))
A=sp.diff(ft(ss,tt,a,b,c)[0],ss)
B=sp.diff(ft(ss,tt,a,b,c)[0],tt)
C=sp.diff(ft(ss,tt,a,b,c)[1],ss)
D=sp.diff(ft(ss,tt,a,b,c)[1],tt)
A*D-B*C
sp.simplify(A*D-B*C)
sp.plot(ss**2)
ft(ss,tt,a,b,c)
sp.simplify(ft(ss,tt,a,b,c))
sp.simplify(ft(ss,tt,a,b,c)[0])
sp.simplify((1,2))
sp.simplify((ss,tt))
sp.simplify((ss,tt**2))
L,fs,fr = sp.var('Ls,fs,fr')
w=L/fr
fa=L*fs/(1+fs)
fb=L/(1+fs)
a=w/2
b=fa
a
b
c=L
a
b
c
L,fs,fr = sp.var('L,fs,fr')
c
L,fs,fr = sp.var('L,fs,fr')
w=L/fr
fa=L*fs/(1+fs)
fb=L/(1+fs)
a=w/2
b=fa
c=L
(a,b,c)
ss
sp.simplify(ft(ss,tt,a,b,c))
sp.simplify(ft(ss,tt,a,b,c))[0]
sp.simplify(ft(ss,tt,a,b,c))[1]
J
J=A*D-B*C
J
A=sp.diff(ft(ss,tt,a,b,c)[0],ss)
B=sp.diff(ft(ss,tt,a,b,c)[0],tt)
C=sp.diff(ft(ss,tt,a,b,c)[1],ss)
D=sp.diff(ft(ss,tt,a,b,c)[1],tt)
J=A*D-B*C
J
sp.pretty
sp.pretty(J)
J.simplify()
sp.pretty(J.simplify())
print(sp.pretty(J.simplify()))
sp.simplify(ft(ss,tt,a,b,c))[1].pretty()
sp.simplify(ft(ss,tt,a,b,c))[1].print()
print(sp.pretty(sp.simplify(ft(ss,tt,a,b,c))))
fr
def rotate(x,y,theta):
    RR = array([[cos(theta),-sin(theta)],[sin(theta),cos(theta)]])
    return RR @ array([[x,y]]).T
rotate(1,0,pi/2)
rotate([1,0],[0,1],pi/2)
x
y
x,y=ft(s,t,2,1,5)
figure(3)
clf();g(*ft(s,t,10,1,5))
clf();g(*ft(x,y,10,1,5))
x,y=ft(s,t,2,1,5)
x
s
t
w
q
s,t=meshgrid(q,q)
x,y=ft(s,t,2,1,5)
clf();g(*ft(x,y,10,1,5))
clf();g(*ft(s,t,10,1,5))
clf();g(rotate(*ft(s,t,2,1,5),pi/2)))
clf();g(rotate(*ft(s,t,2,1,5),pi/2))
clf();g(*rotate(*ft(s,t,2,1,5),pi/2))
x
rotate(x,y,pi/2)
plot(rotate(x,y,pi/2),'o-
plot(rotate(x,y,pi/2),'o-')
xp,yp=rotate(x,y,pi/2)
aa=rotate(x,y,pi/2)
aa
aa.shape
def rotate(x,y,theta):
    RR = array([[cos(theta),-sin(theta)],[sin(theta),cos(theta)]])
    return RR @ array([[x,y]]).T
def ft(s,t,a,b,c,theta):
    x = a/4*((s+1)*(1-t)+(t+1)*(s-1))
    y = b/4*((t+1)*(1-s)+(s+1)*(2+(t+1)*(c/b-1)))
    return rotate(x,y,theta-pi/2)
ft(ss,tt,a,b,c,0)
sp.Matrix
def srotate(x,y,theta):
    RR = sp.Matrix([[sp.cos(theta),-sp.sin(theta)],[sp.sin(theta),sp.cos(theta)]])
    return RR * sp.Matrix([[x,y]]).T
sp.Matrix([[1],[0]])
def ft(s,t,a,b,c,theta):
    x = a/4*((s+1)*(1-t)+(t+1)*(s-1))
    y = b/4*((t+1)*(1-s)+(s+1)*(2+(t+1)*(c/b-1)))
    return srotate(x,y,theta-pi/2)
ft(ss,tt,a,b,c,0)
ft(ss,tt,a,b,c,tt)
tt
theta
theta=sp.var('theta')
ft(ss,tt,a,b,c,theta)
def ft(s,t,a,b,c,theta):
    x = a/4*((s+1)*(1-t)+(t+1)*(s-1))
    y = b/4*((t+1)*(1-s)+(s+1)*(2+(t+1)*(c/b-1)))
    return srotate(x,y,theta-sp.pi/2)
ft(ss,tt,a,b,c,theta)
sp.simplify(ft(ss,tt,a,b,c,theta))
print(sp.pretty(sp.simplify(ft(ss,tt,a,b,c,theta))))
A=sp.diff(ft(ss,tt,a,b,c)[0],ss)
A=sp.diff(ft(ss,tt,a,b,c,theta)[0],ss)
B=sp.diff(ft(ss,tt,a,b,c,theta)[0],tt)
C=sp.diff(ft(ss,tt,a,b,c,theta)[1],ss)
D=sp.diff(ft(ss,tt,a,b,c,theta)[1],tt)
J=A*D-B*C
J=(A*D-B*C).simplify()
J
print(sp.pretty(J))
get_ipython().magic('save frond_transform_algebra')
get_ipython().magic('save frond_transform_algebra *')
get_ipython().magic('save frond_transform_algebra 0-224')
