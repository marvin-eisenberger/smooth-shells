function Tblock = translationtoblocks(Tw,sizeM)

    Tblock = kron(spdiags(Tw(1,:)',0,sizeM,sizeM),[zeros(4,3),[1;0;0;0]]) + kron(spdiags(Tw(2,:)',0,sizeM,sizeM),[zeros(4,3),[0;1;0;0]]) + kron(spdiags(Tw(3,:)',0,sizeM,sizeM),[zeros(4,3),[0;0;1;0]]);

end

