import numpy as np 
### Implement the Fourier-Cosine Method
def HestonCOS(S,K,T,r,sigma,lmbda,meanV,v0,rho, otype, N=256):
    c1 = r*T+(1-np.exp(-lmbda*T))*(meanV-v0)/(2.0*lmbda)-0.5*meanV*T
    c2 = 1.0/(8.0*lmbda**3)*(sigma*T*lmbda*np.exp(-lmbda*T)*(v0-meanV)*(8.0*lmbda*rho-4.0*sigma)+lmbda*rho*sigma*(1-np.exp(-lmbda*T))*(16.0*meanV-8.0*v0)+2.0*meanV*lmbda*T*(-4.0*lmbda*rho*sigma+sigma**2+4.0*lmbda**2)+sigma**2*((meanV-2.0*v0)*np.exp(-2.0*lmbda*T)+meanV*(6.0*np.exp(-lmbda*T)-7.0)+2.0*v0)+8.0*lmbda**2*(v0-meanV)*(1-np.exp(-lmbda*T)))
    a = c1-12.0*np.sqrt(np.abs(c2))
    b = c1+12.0*np.sqrt(np.abs(c2))
    x = np.log(np.array(S, dtype='float')/K)
    k = np.arange(0,N)
    if (otype == 'C'):
       U = 2.0/(b-a)*(xi(k,a,b,0,b) - psi(k,a,b,0,b))
    else:
       U = 2.0/(b-a)*(-xi(k,a,b,a,0) + psi(k,a,b,a,0))   
    unit = [1.0] * N
    unit[0] = 0.5
    ret = 0
# Note that HestonCF is independent of the strike
    HCF = HestonCF(k*np.pi/(b-a),T,r,sigma,lmbda,meanV,v0,rho)
    for i in range(N):
      ret += unit[i]*HCF[i]*np.exp(1j*np.float(k[i])*np.pi*(x-a)/(b-a))*U[i]
    return K*np.exp(-r*T)*ret.real

def HestonVega(S,K,T,r,sigma,lmbda,meanV,v0,rho, otype, N=256):
    c1 = r*T+(1-np.exp(-lmbda*T))*(meanV-v0)/(2.0*lmbda)-0.5*meanV*T
    c2 = 1.0/(8.0*lmbda**3)*(sigma*T*lmbda*np.exp(-lmbda*T)*(v0-meanV)*(8.0*lmbda*rho-4.0*sigma)+lmbda*rho*sigma*(1-np.exp(-lmbda*T))*(16.0*meanV-8.0*v0)+2.0*meanV*lmbda*T*(-4.0*lmbda*rho*sigma+sigma**2+4.0*lmbda**2)+sigma**2*((meanV-2.0*v0)*np.exp(-2.0*lmbda*T)+meanV*(6.0*np.exp(-lmbda*T)-7.0)+2.0*v0)+8.0*lmbda**2*(v0-meanV)*(1-np.exp(-lmbda*T)))
    a = c1-12.0*np.sqrt(abs(c2))
    b = c1+12.0*np.sqrt(abs(c2))
    x = np.log(np.array(S, dtype="float")/K)
    k = np.arange(0,N)
    if (otype == 'C'):
       U = 2.0/(b-a)*(xi(k,a,b,0,b) - psi(k,a,b,0,b))
    else:
       U = 2.0/(b-a)*(-xi(k,a,b,a,0) + psi(k,a,b,a,0))   
    unit = [1.0] * N
    unit[0] = 0.5
    ret = 0
# Note that HestonCF is independent of the strike
    HCF = HestonCFdu0(k*np.pi/(b-a),T,r,sigma,lmbda,meanV,v0,rho)
    for i in range(N):
      ret += unit[i]*HCF[i]*np.exp(1j*np.float(k[i])*np.pi*(x-a)/(b-a))*U[i]
    return K*np.exp(-r*T)*ret.real


def HestonCF(u,T,r,sigma,lmbda,meanV,v0,rho):
    a = lmbda*meanV
    b = lmbda
    d = np.sqrt((1j*rho*sigma*u-b)**2+(u**2+1j*u)*sigma**2)
    g = (b-1j*rho*sigma*u-d)/(b-1j*rho*sigma*u+d)
    ret = np.exp(1j*u*r*T)
    ret = ret*np.exp((a/sigma**2)*((b - rho*1j*sigma*u - d)*T - 2.0*np.log((1-g*np.exp(-d*T))/(1-g))))
    return ret*np.exp((v0/sigma**2)*(b - rho*1j*sigma*u - d)*(1-np.exp(-d*T))/(1-g*np.exp(-d*T)))


# The Derivative of the Heston Characteristic Function w.r.t. to the v0 
def HestonCFdu0(u,T,r,sigma,lmbda,meanV,v0,rho):

    a = lmbda*meanV
    b = lmbda
    d = np.sqrt((1j*rho*sigma*u-b)**2+(u**2+1j*u)*sigma**2)
    g = (b-1j*rho*sigma*u-d)/(b-1j*rho*sigma*u+d)
    ret = np.exp(1j*u*r*T)
    
    ret = ret*np.exp((a/sigma**2)*((b - rho*1j*sigma*u - d)*T - 2.0*np.log((1.0-g*np.exp(-d*T))/(1.0-g))))
    ret = ret*np.exp((v0/sigma**2)*(b - rho*1j*sigma*u - d)*(1.0-np.exp(-d*T))/(1.0-g*np.exp(-d*T)))
    ret = ret *((1.0-np.exp(-d*T))/(1.0-g*np.exp(-d*T)))*(b-1j*rho*sigma*u-d)/(sigma**2)
    return(ret)



def xi(k,a,b,c,d):
    ret = 1.0/(1+(k*np.pi/(b-a))**2)*(np.cos(k*np.pi*(d-a)/(b-a))*np.exp(d)-np.cos(k*np.pi*(c-a)/(b-a))*np.exp(c)+k*np.pi/(b-a)*np.sin(k*np.pi*(d-a)/(b-a))*np.exp(d)-k*np.pi/(b-a)*np.sin(k*np.pi*(c-a)/(b-a))*np.exp(c))
    return ret

def psi(k,a,b,c,d):
    N = len(k)
    idx = np.arange(1, N)
    ret = np.array([0.0]*N)
    ret[0] = d-c
    ret[idx] =(np.sin(k[idx]*np.pi*(d-a)/(b-a))-np.sin(k[idx]*np.pi*(c-a)/(b-a)))*(b-a)/(k[idx]*np.pi)
    return ret
