% Jacobian of explicit kinematic constraints of
% fourbar2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
% 
% Output:
% W [5x1]
%  Derivative of the joint coordinates w.r.t minimal coordinates
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:27
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function W = fourbar2DE2_kinconstr_expl_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar2DE2_kinconstr_expl_jacobian_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2DE2_kinconstr_expl_jacobian_mdh_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:27:44
% EndTime: 2020-04-24 20:27:44
% DurationCPUTime: 0.01s
% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [1; -1; 1; 1; 0;];
W = t1;
