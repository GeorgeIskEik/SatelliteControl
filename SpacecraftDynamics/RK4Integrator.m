clear all 
% parameter initialization
y0 = [0.408248;0.;0.408248;0.816497];
t0 = 1;
tf = 60;
dx = 0.01;
t = t0:dx:tf;
yn = y0;
store = zeros(4,length(t));

for i = 1:length(t)
    % calculate the 4 slope weighed increments.
    k1 = Omega(i,yn(:,i));
    k2 = Omega(i+1/2*dx,yn(:,i)+1/2*dx*k1);
    k3 = Omega(i+1/2*dx,yn(:,i)+1/2*dx*k2);
    k4 = Omega(i+dx,yn(:,i)+dx*k3);
    % update vector value.
    yn(:,i+1) = yn(:,i)+dx/6*(k1+k2+k3+k4);
    % update external storage.
    store(:,i) = yn(:,i);
end
% pul out the quaternion parameters
b0 = store(1,:);
b1 = store(2,:);
b2 = store(3,:);
b3 = store(4,:);

% plot the quaternion parameters
figure 
subplot(4,1,1)
plot(t,(2*acos(b0)),'linewidth',2)
subplot(4,1,2)
plot(t,b1,'linewidth',2)
subplot(4,1,3)
plot(t,b2,'linewidth',2)
subplot(4,1,4)
plot(t,b3,'linewidth',2)

% extract the norm of the vectorial part for a sample @ t = 42(s)
% Find the index corresponding to the closest time value to 42 seconds
[~, sample_index] = min(abs(t - 42));

% Calculate the norm of the vectorial part at the found index
norm = sqrt(store(2,sample_index).^2 + store(3,sample_index).^2 + store(4,sample_index).^2);

% Display the calculated norm
disp(["The norm @ 42s is: ", num2str(norm)])


% sample_index = find(abs(t - 42) < 1e-10, 1);
% norm = sqrt(store(2,sample_index).^2+store(3,sample_index).^2+store(4,sample_index).^2);
% disp(["The norm @ 42s is: ",num2str(norm)])