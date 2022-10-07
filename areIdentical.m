function result = areIdentical(lhs, rhs)
  if size(lhs) ~= size(rhs)
    result = false;
  else
    result = isequal(isAlways(lhs == rhs), ones(size(lhs)));
  end
end
