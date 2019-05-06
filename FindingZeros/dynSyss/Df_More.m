function J = Df_More(x)
J = zeros(length(x),length(x));
% for i = 1:size(J,1)
%     for j = 1:size(J,2)
%         if i == j
%             J(i,j) = -sin(x(j))+1-cos(x(i))+i*sin(x(i))-cos(x(i));
%         else
%             J(i,j) = -sin(x(j));
%         end
%     end
% end
%
% squeeze every possible chance to vectorize my code
J = repmat(-sin(x),1,length(x));
J = J + diag( 1-2*cos(x)+diag(1:length(x))*sin(x) );