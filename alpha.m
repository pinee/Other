function f=alpha(u1,u2,E1,E2)
if u2=u1
    f=avg(u1,u2);
else
    f=(E2-E1)/(u2-u1);
end