function p = penaltymin(val,minval,dist)

p = max(tanh((minval-val)./(0.5*dist)),0);

