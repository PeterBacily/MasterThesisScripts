from __future__ import division


def poker(pot_flop,betsize1,betsize2,betsize3):
    tp = (1+2*betsize1)*pot_flop
    print 'tp=', tp, betsize1*pot_flop
    rp = (1+2*betsize2)*tp
    print 'rp =',rp,betsize2*tp
    fp =(1+2*betsize3)*rp
    print 'fp = ',fp,betsize3*rp

poker(7.5,2/3,2/3,2/3)

print '-----'

poker(7.5,3/4,3/4,3/4)

print '_____'

poker(7.5,1,0.8,3/4)