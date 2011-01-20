% vertex transition trajectories

tmax = 46;

t1 = [4 4 4 5 4 5 5 5 5 5 4 5 5 5 6 7 7 8 8 9 9 9 9 10];
length(t1)
t2 = [2 2 2 2 2 2 2 2 3 2 2 2 2 2 2 3 3 3 4 2 3 3 3 2 2 2 3 3 2 3 5 5 5 5 5 8 10 10 10 11 11];
length(t2)
t3 = [4 4 4 5 3 5 5 5 5 5 5 4 4 5 5 5 5 5 5 6 6 6 6 8 8 9];
length(t3)
t4 = [6 6 6 6 6 6 6 6 7 7 7 7 9 9 8 10 10 10 10 10];
length(t4)
t8 = [3 3 3 3 3 3 4 4 5 4 4 4 5 4 4 4 4 4 5 5 5 5 5 4 4 5 5 4 5 5 5 5 6 5 5 6 6 6 7 8 8 8 9 9];
length(t8)


t1 = pad_with_NaN(t1, tmax);
t2 = pad_with_NaN(t2, tmax);
t3 = pad_with_NaN(t3, tmax);
t4 = pad_with_NaN(t4, tmax);
t8 = pad_with_NaN(t8, tmax);


t = [t1 t2 t3 t4 t8];

figure;
plot((0:45)*19.5/60,-2.5*t)
xlabel('time (min)');
ylabel('distance from uppermost plane (\mu m)');
