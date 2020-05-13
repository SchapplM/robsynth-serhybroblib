% Jacobian of explicit kinematic constraints of
% fourbarprisDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% W [5x1]
%  Derivative of the joint coordinates w.r.t minimal coordinates
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:10
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function W = fourbarprisDE1_kinconstr_expl_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_kinconstr_expl_jacobian_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_kinconstr_expl_jacobian_mdh_sym_varpar: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:08:36
% EndTime: 2020-05-07 09:08:36
% DurationCPUTime: 0.02s
% Computational Cost: add. (50->13), mult. (28->13), div. (2->1), fcn. (6->2), ass. (0->8)
t28 = (qJ(1) + pkin(3));
t34 = -pkin(2) + t28;
t35 = -pkin(2) - t28;
t26 = ((-(pkin(1) + t35) * (pkin(1) + t34) * (pkin(1) - t34) * (pkin(1) - t35)) ^ (-0.1e1 / 0.2e1));
t37 = (t26 / t28);
t36 = (pkin(1) ^ 2 - pkin(2) ^ 2);
t33 = -qJ(1) ^ 2 + (-2 * qJ(1) - pkin(3)) * pkin(3);
t1 = [(t33 + t36) * t37; 1; -2 * t28 * t26; (-t33 + t36) * t37; 0;];
W = t1;
