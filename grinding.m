
iGood = 0;

while true
  example_script;
  if mean(normv(vertCurrFull-Y.vert)) > 0.025
    break;
  end
  iGood = iGood+1
end
