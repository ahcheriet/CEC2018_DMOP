# The code is translated from Matlab by ahcheriet@gmail.com
# ========================================================|#
# The 14 test functions are for cec2018 competition on    |
# dynamic multiobjective optimisation. This document is   |
# free to disseminate for academic use.                   |
# --------------------------------------------------------|#
# The "time" term in the test suite is defined as:        |
#          t=1/nt*floor(tau/taut)                         |
# where - nt:    severity of change                       |
#       - taut:  frequency of change                      |
#       - tau:   current generation counter               |
# --------------------------------------------------------|#
# Any questions can be directed to                        |
#    Dr. Shouyong Jiang at math4neu@gmail.com.            |
#                                                         |
# ========================================================|#

# cec2018_DF(probID, x, tau, taut, nt)
# INPUT:
#       probID: test problem identifier (i.e. 'DF1')
#       x:      variable vector
#       tau:    current generation counter
#       taut:   frequency of change
#       nt:     severity of change
#
# OUTPUT:
#       f:      objective vector
#
from numpy import sin,cos,power,setdiff1d,exp,floor,pi,arange,prod



def cec2018_DF(problemID=None,x=None,tau=None,taut=None,nt=None):
# INPUT:
#       probID: test problem identifier (i.e. 'DF1')
#       x:      variable vector
#       tau:    current generation counter
#       taut:   frequency of change
#       nt:     severity of change
    
# OUTPUT:
#       f:      objective vector
    
# the first change occurs after T0 generations, that is, the
# generation at which a change occurs is (T0+1), (T0+taut+1), etc.

    T0=50
    # calculate time instant
    tau_tmp=max(tau + taut - (T0 + 1),0)
    t=(1.0 / nt)*floor(tau_tmp / taut)
    n=len(x)
    f={}
    if problemID == 'DF1':
        G=abs(sin(0.5*pi*t))
        H=0.75*sin(0.5*pi*t)+1.25
        g=1+sum((x[1:]-G)**2)
        f[0]=x[0]
        f[1]=g*power(1-(x[0]/g),H)
    if problemID == 'DF2':
        G=abs(sin(0.5*pi*t))
        r=1+floor((n-1)*G)
        tmp=setdiff1d(range(0,n),[int(r)])
        g=1+sum([(x[int(index)]-G)**2 for index in tmp])
        f[0]=x[int(r)]
        f[1]=g*(power(1-(x[int(r)]/g),0.5))
    if problemID == 'DF3':
        G=sin(0.5*pi*t)
        H=G+1.5
        g=1+sum( power(x[1:]-G-x[0],H)**2)
        f[0]=x[0]
        f[1]=g*power(1-(x[0]/g),H)
    if problemID == 'DF4':
        a=sin(0.5*pi*t)
        b=1+abs(cos(0.5*pi*t))
        H=1.5+a
        g=1+sum((x[1:]-a*x[0]**2/x[1:])**2)
        f[0]=g*power(abs(x[0]-a),H)
        f[1]=g*power(abs(x[0]-a-b),H)
    if problemID == 'DF5':
        G=sin(0.5*pi*t)
        w=floor(10*G)
        g=1+sum((x[1:]-G)**2)
        f[0]=g*(x[0]+0.02*sin(w*pi*x[0]))
        f[1]=g*(1-x[0]+0.02*sin(w*pi*x[0]))
    if problemID == 'DF6':
        G=sin(0.5*pi*t)
        a=0.2+2.8*abs(G)
        y=x[1:]-G
        g=1+sum((abs(G)*y**2-10*cos(2*pi*y)+10))
        f[0]=g*power(x[0]+0.1*sin(3*pi*x[0]),a)
        f[1]=g*power(1-x[0]+0.1*sin(3*pi*x[0]),a)
    if problemID == 'DF7':
        a=5*cos(0.5*pi*t)
        tmp=1/(1+exp(a*(x[0]-2.5)))
        g=1+sum(power(x[1:]-tmp,2))
        f[0]=g*(1+t)/x[0]
        f[1]=g*x[0]/(1+t)        
    if problemID == 'DF8':
        G=sin(0.5*pi*t)
        a=2.25+2*cos(2*pi*t)
        b=100*G**2
        tmp=G*sin(power(4*pi*x[0],b))/(1+abs(G))
        g=1+sum((x[1:]-tmp)**2)
        f[0]=g*(x[0]+0.1*sin(3*pi*x[0]))
        f[1]=g*power(1-x[1]+0.1*sin(3*pi*x[1]),a)        
    if problemID == 'DF9':
        N=1+floor(10*abs(sin(0.5*pi*t)))
        g=1
        for i in range(1,n):
            tmp=x[i]-cos(4*t+x[0]+x[i-1])
            g=g+tmp**2
        f[0]=g*(x[0]+max(0, (0.1+0.5/N)*sin(2*N*pi*x[0])))
        f[1]=g*(1-x[0]+max(0, (0.1+0.5/N)*sin(2*N*pi*x[0])))        
    if problemID == 'DF10':
        G=sin(0.5*pi*t)
        H=2.25+2*cos(0.5*pi*t)
        tmp=sin(2*pi*(x[0]+x[1]))/(1+abs(G))
        g=1+sum((x[2:]-tmp)**2)
        f[0]=g*power(sin(0.5*pi*x[0]),H)
        f[1]=g*power(sin(0.5*pi*x[1]),H)*power(cos(0.5*pi*x[0]),H)
        f[2]=g*power(cos(0.5*pi*x[1]),H)*power(cos(0.5*pi*x[0]),H)
    if problemID == 'DF11':
        G=abs(sin(0.5*pi*t))
        g=1+G+sum((x[2:]-0.5*G*x[0])**2)
        y=[pi*G/6.0+(pi/2-pi*G/3.0)*x[i] for i in [0,1]]
        f[0]=g*sin(y[0]) 
        f[1]=g*sin(y[1])*cos(y[0])
        f[2]=g*cos(y[1])*cos(y[0])
    if problemID == 'DF12':
        k=10*sin(pi*t)
        tmp1=x[2:]-sin(t*x[0])
        tmp2=[sin(floor(k*(2*x[0]-1))*pi/2)]
        g=1+sum(tmp1**2)+prod(tmp2)
        f[0]=g*cos(0.5*pi*x[1])*cos(0.5*pi*x[0])
        f[1]=g*sin(0.5*pi*x[1])*cos(0.5*pi*x[0])
        f[2]=g*sin(0.5*pi*x[1])
    if problemID == 'DF13':
        G=sin(0.5*pi*t);
        p=floor(6*G);
        g=1+sum((x[2:]-G)**2)
        f[0]=g*cos(0.5*pi*x[0])**2
        f[1]=g*cos(0.5*pi*x[1])**2
        f[2]=g*sin(0.5*pi*x[0])**2+sin(0.5*pi*x[0])*cos(p*pi*x[0])**2+sin(0.5*pi*x[1])**2+sin(0.5*pi*x[1])*cos(p*pi*x[1])**2
    if problemID == 'DF14':
        G=sin(0.5*pi*t)
        g=1+sum((x[2:]-G)**2)
        y=0.5+G*(x[0]-0.5)
        f[0]=g*(1-y+0.05*sin(6*pi*y))
        f[1]=g*(1-x[1]+0.05*sin(6*pi*x[1]))*(y+0.05*sin(6*pi*y))
        f[2]=g*(x[1]+0.05*sin(6*pi*x[1]))*(y+0.05*sin(6*pi*y))
    return f
if __name__ == '__main__':
    for i in range(1,15):
        pb = 'DF'+str(i)
        print pb
        print cec2018_DF(problemID=pb,x=[0.1,0.2,0.3,0.4],tau=1,taut=1,nt=1)
#        except Exception,e:
#            print str(e)
#            print '============='
#            pass
    

