import numpy as np
import math
from math import sin, cos, tan, atan
import matplotlib.pyplot as plt


#############################################################################

a = 120 # mm Abstand Rad Loslager A
b = 170 # mm Abstand Rad Festlager B
c = 60 # mm Abstand Rad Loslager C
e = 150 # mm Abstand Loslager C Festlager E
d = 80 # mm Durchmesser Tellerrad
r = d/2 # mm Radius Tellerrad
i = 2 # Uebersetzung
z1 = 50 # Zähneanzahl Tellerrad
z2 = z1/i # Zähneanzahl Ritzel
achsenwinkel = 90 # ° Achswinkel
betta = 20 # ° Schraegungswinkel
alpha = 20 # ° Eingriffswinkel
p = 15000 # W Leistun an der Welle
n = 1500 # 1/min Umdrehungen
k_a = 1.5 # Betriebsfaktor nach der Tabelle
alpha_0 = 0.7 # Anstrengungsverhaeltnis nach der Tabelle
sigma_zul = 245 # MPa Zulaessige Spannung nach der Tabelle
r2 = r/i # Radius Ritzel

## Einige Hilfsrechnungen /Umrechnungen
achsenwinkel = achsenwinkel*math.pi/180
u = z1/z2
delta1 = atan(sin(achsenwinkel)/(u + cos(achsenwinkel)))
#delta1*180/math.pi
delta2 = achsenwinkel - delta1
betta = betta*math.pi/180
alpha = alpha*math.pi/180
n = n/60
t_moment = p*k_a/(2*math.pi*n)
f_t = t_moment/r*1000
a = a/1000 # m
b = b/1000 # m
c = c/1000 # m
e = e/1000 # m
r = r/1000 # m

#############################################################################

f_a1 = (sin(delta1)*tan(alpha)+cos(delta1)*sin(betta))*f_t/cos(betta)
f_a2 = (sin(delta2)*tan(alpha)-cos(delta2)*sin(betta))*f_t/cos(betta)
f_r1 = (cos(delta1)*tan(alpha)-sin(delta1)*sin(betta))*f_t/cos(betta)
f_r2 = (cos(delta2)*tan(alpha)+sin(delta2)*sin(betta))*f_t/cos(betta)
f_t1 = f_t
f_t2 = f_t
print('Tangentiale Kraft Tellerrad =', round(f_t1,1), 'N') #
print('Radiale Kraft Tellerrad =', round(f_r1,1), 'N')
print('Axiale Kraft Tellerrad =', round(f_a1,1), 'N')
print('Tangentiale Kraft Ritzel =', round(f_t2,1), 'N')
print('Radiale Kraft Ritzel =', round(f_r2,1), 'N')
print('Axiale Kraft Ritzel =', round(f_a2,1), 'N')


#############################################################################

A = [[0,0,-1,0,0], [1,0,0,1,0],[0,-1,0,0,-1], [0,-(a+b), 0, 0, 0], [-(a+b), 0,0,0,0]]
vb = [-f_a1,f_r1, -f_t1, -b*f_t1, -(b*f_r1-r*f_a1)]

#############################################################################

res = np.linalg.solve(A, vb)
print('Y Kraft Lager A =', round(res[0],1), 'N')
print('Z Kraft Lager A =', round(res[1],1), 'N')
print('X Kraft Lager B =', round(res[2],1), 'N')
print('Y Kraft Lager B =', round(res[3],1), 'N')
print('Z Kraft Lager B =', round(res[4],1), 'N')
#np.allclose(np.dot(A, res), vb)
x, m1, m2, mg = [], [],[], []
i = 0.01
f_a = res[0]
while i <= (a+b):
    x.append(i)
    i+=0.001
for xp in x:
    if xp < a:
        mq = res[0]*xp
        m1.append(mq)
        mn = res[1]*xp
        m2.append(mn)
        mg.append(math.sqrt(mn**2 + mq**2))
    else:
        #m1 = -res[3]*(b-xp+a)
        mq = -f_r1*(xp-a)+f_a1*r+res[0]*xp
        m1.append(mq)
        mn = res[1]*xp - f_t1*(xp-a)
        m2.append(mn)
        mg.append(math.sqrt(mn**2 + mq**2))
        
        #m2 = -res[4]*(b-(xp+a))
plt.plot(x, m1, color = 'orange')
plt.plot(x, m2, color = 'red')
plt.plot(x, mg, color = 'blue')
plt.ylabel('Biegemoment [Nm]')
plt.xlabel('Abstand vom Lager A')
plt.vlines(0, ymin = 0, ymax = 5)
plt.vlines(a+b, ymin = 0, ymax = 5)
plt.hlines(0, xmin = 0, xmax = a+b)
plt.show()

#############################################################################

f_a_radial = math.sqrt(res[0]**2 + res[1]**2)
f_b_radial = math.sqrt(res[3]**2 + res[4]**2)
print('Radiale Kraft Lager A =', round(f_a_radial,1), 'N')
print('Radiale Kraft Lager B =', round(f_b_radial,1), 'N')


#############################################################################

m_max_z = abs(-res[0]*a)

if m_max_z < abs(res[3]*b):
    m_max_z = abs(res[3]*b)

m_max_y = abs(-res[1]*a)
if m_max_y < abs(res[4]*b):
    m_max_y = abs(res[4]*b)
print ('Maximales Biegemoment in xz-Ebene:', round(m_max_y,1),'Nm')
print ('Maximales Biegemoment in xy-Ebene:', round(m_max_z,1),'Nm')

############################################################################

b_moment= math.sqrt(m_max_y**2 + m_max_z**2)
print('Maximales Biegemoment =', round(b_moment,1), 'Nm')

print(-res[0]*a,-res[3]*b)
print(-res[1]*a,-res[4]*b)
print(f_a1*r)

v_moment=math.sqrt(b_moment**2 + 0.75*(alpha_0*t_moment)**2)
print('Vergleichsmoment =', round(v_moment,1),'Nm')

d_erf=(32*v_moment*1000/math.pi/sigma_zul)**(1/3)
print('erforderlicher Durchmesser =', round(d_erf,1), 'mm')


############################################################################
############################################################################
############################################################################
############################################################################

# Berechnung der Ritzelwelle
r2 = r2/1000
###########################
f_a1 = (sin(delta1)*tan(alpha)+cos(delta1)*sin(betta))*f_t/cos(betta)
f_a2 = (sin(delta2)*tan(alpha)-cos(delta2)*sin(betta))*f_t/cos(betta)
f_r1 = (cos(delta1)*tan(alpha)-sin(delta1)*sin(betta))*f_t/cos(betta)
f_r2 = (cos(delta2)*tan(alpha)+sin(delta2)*sin(betta))*f_t/cos(betta)
f_t1 = f_t
f_t2 = f_t
print('Tangentiale Kraft Ritzel =', round(f_t2,1), 'N')
print('Radiale Kraft Ritzel =', round(f_r2,1), 'N')
print('Axiale Kraft Ritzel =', round(f_a2,1), 'N')

A = [[-1,0,-1,0,0], [0,0,0,-1,0], [0,-1,0,0,-1], [0,0,e,0,0], [0,0,0,0,e]]
vb = [-f_r2, -f_a2, -f_t2, -c*f_r2+r2*f_a2, -c*f_t2]
res = np.linalg.solve(A, vb)
print('X Kraft Lager A =', round(res[0],1), 'N')
print('Z Kraft Lager A =', round(res[1],1), 'N')
print('X Kraft Lager B =', round(res[2],1), 'N')
print('Y Kraft Lager B =', round(res[3],1), 'N')
print('Z Kraft Lager B =', round(res[4],1), 'N')
print('##############################################################')


x, m1, m2, mg = [], [],[], []
i = 0.0
while i <= (c+e):
    x.append(i)
    i+=0.001
for xp in x:
    if xp <= e:
        mq=(-res[2]*(xp))
        m1.append(mq)
        mn=(-res[4]*(xp))
        m2.append(mn)
        mg.append(math.sqrt(mn**2 + mq**2))
    else:
        #m1.append(-f_r2*(xp-a)+f_a2*r+res[0]*xp)
        #m2.append(res[1]*xp - f_t2*(xp-a))
        
        mq=(-res[2]*(xp)-res[0]*(xp-e))
        m1.append(mq)
        mn=-res[4]*(xp)-res[1]*(xp-e)
        m2.append(mn)
        mg.append(math.sqrt(mn**2 + mq**2))
plt.plot(x, m1, color = 'orange')
plt.plot(x, m2, color = 'red')
plt.plot(x, mg, color = 'blue')
plt.ylabel('Biegemoment [Nm]')
plt.xlabel('Abstand vom Lager E')
plt.vlines(0, ymin = 0, ymax = 5)
plt.vlines(e+c, ymin = 0, ymax = 5)
plt.hlines(0, xmin = 0, xmax = e+c)
plt.show()

print('##############################################################')


f_a_radial = math.sqrt(res[0]**2+res[1]**2)
f_b_radial = math.sqrt(res[2]**2+res[4]**2)
print('Radiale Kraft Lager A =', round(f_a_radial,1), 'N')
print('Radiale Kraft Lager B =', round(f_b_radial,1), 'N')
print('##############################################################')
m_max_z = abs(f_r2*c)
if m_max_z < abs(res[2]*e):
    m_max_z = abs(res[2]*e)
m_max_y = abs(f_a2*c)
if m_max_y < abs(res[4]*e):
    m_max_y = abs(res[4]*e)
print ('Maximales Biegemoment in xz-Ebene:', round(m_max_y,1),'Nm')
print ('Maximales Biegemoment in xy-Ebene:', round(m_max_z,1),'Nm')
print('##############################################################')
b_moment= math.sqrt(m_max_y**2 + m_max_z**2)
print('Maximales Biegemoment =', round(b_moment,1), 'Nm')
print('##############################################################')
print(f_r2*c,res[2]*e)
print(f_t2*c,res[4]*e)
print(f_a1*r)
print('##############################################################')
v_moment=math.sqrt(b_moment**2 + 0.75*(alpha_0*t_moment)**2)
print('Vergleichsmoment =', round(v_moment,1),'Nm')
print('##############################################################')
#sigma_zul = sigma_zul/4
d_erf=(32*v_moment*1000/math.pi/sigma_zul)**(1/3)
print('erforderlicher Durchmesser =', round(d_erf,1), 'mm')
print ('Mittlerer Kegelraddurchmesser:', r2*2*1000, 'mm')
