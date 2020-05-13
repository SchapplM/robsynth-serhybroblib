% Jacobian of explicit kinematic constraints of
% fourbar1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% W [5x1]
%  Derivative of the joint coordinates w.r.t minimal coordinates
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:49
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function W = fourbar1TE_kinconstr_expl_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_kinconstr_expl_jacobian_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_kinconstr_expl_jacobian_mdh_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:46:44
% EndTime: 2020-04-24 19:46:44
% DurationCPUTime: 0.09s
% Computational Cost: add. (82->19), mult. (100->24), div. (2->2), fcn. (29->4), ass. (0->13)
t48 = -pkin(3) - pkin(4);
t47 = -pkin(3) + pkin(4);
t34 = sin(qJ(1));
t46 = pkin(2) * t34;
t35 = cos(qJ(1));
t45 = pkin(2) * t35;
t43 = (-0.2e1 * t45 + pkin(1)) * pkin(1);
t40 = sqrt(-((pkin(2) - t48) * (pkin(2) + t48) + t43) * ((pkin(2) - t47) * (pkin(2) + t47) + t43));
t29 = 0.1e1 / t40;
t41 = pkin(2) ^ 2 + t43;
t44 = t29 / t41;
t42 = pkin(3) ^ 2 - pkin(4) ^ 2;
t1 = [1; -pkin(1) * ((pkin(1) - t45) * t40 + (t41 - t42) * t46) * t44; pkin(2) * ((-pkin(1) * t35 + pkin(2)) * t40 + pkin(1) * t34 * (t41 + t42)) * t44; 0.2e1 * pkin(1) * t29 * t46; 0;];
W = t1;
