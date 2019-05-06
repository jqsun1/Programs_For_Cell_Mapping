function J = Df_neural(f,X)
dX = 0.001*ones(length(X),1);
X_aug = repmat(X',length(X),1) + diag(dX);
fx = f(X);
J = [(f(X_aug(1,:)') - fx)/dX(1),...
     (f(X_aug(2,:)') - fx)/dX(2)];