% Jacobian of explicit kinematic constraints of
% fourbar1turnDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% W [6x2]
%  Derivative of the joint coordinates w.r.t minimal coordinates
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function W = fourbar1turnDE1_kinconstr_expl_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_kinconstr_expl_jacobian_mdh_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_kinconstr_expl_jacobian_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:25:18
% EndTime: 2020-04-12 19:25:18
% DurationCPUTime: 0.28s
% Computational Cost: add. (646->46), mult. (865->83), div. (34->12), fcn. (233->4), ass. (0->46)
t82 = -2 * pkin(2);
t81 = (-pkin(3) - pkin(4));
t80 = (-pkin(3) + pkin(4));
t55 = sin(qJ(2));
t79 = pkin(2) * t55;
t78 = t55 * pkin(1);
t56 = cos(qJ(2));
t77 = t56 * pkin(1);
t69 = pkin(2) * t77;
t53 = -0.2e1 * t69;
t62 = pkin(1) ^ 2;
t72 = t53 + t62;
t42 = ((pkin(2) - t81) * (pkin(2) + t81)) + t72;
t43 = ((pkin(2) - t80) * (pkin(2) + t80)) + t72;
t75 = t42 * t43;
t63 = sqrt(-t75);
t70 = pkin(2) * t78;
t76 = 0.1e1 / t63 * (-t42 - t43) * t70;
t74 = t55 * t63;
t73 = t56 * t63;
t61 = pkin(2) ^ 2;
t71 = t61 + t62;
t68 = pkin(2) * t74;
t50 = t53 + t71;
t49 = 0.1e1 / t50 ^ 2;
t67 = t49 * t79;
t59 = pkin(3) ^ 2;
t66 = -t59 + t71;
t57 = pkin(4) ^ 2;
t54 = t55 ^ 2;
t52 = -pkin(2) + t77;
t51 = -pkin(2) * t56 + pkin(1);
t48 = 0.1e1 / t50;
t47 = t53 + t57 + t66;
t46 = -t57 + t59 + t50;
t45 = t57 - t66 + 0.2e1 * t69;
t44 = 0.1e1 / t45 ^ 2;
t41 = t46 * t78;
t40 = t47 * t79;
t35 = -t52 * t63 + t41;
t34 = t51 * t63 + t40;
t33 = t47 * t51 - t68;
t32 = -pkin(1) * t74 - t46 * t52;
t31 = 0.1e1 / t33 ^ 2;
t30 = 0.1e1 / t32 ^ 2;
t1 = [1, 0; 0, 1; 0, (((0.2e1 * t62 * t54 * pkin(2) - t52 * t76) * t48 + ((t56 * t46 + t74) * t48 - 0.2e1 * t35 * t67) * pkin(1)) / t32 - (t41 * t48 + (-t48 * t73 + ((t52 * t82 - t76) * t48 + t32 * t49 * t82) * t55) * pkin(1)) * t35 * t30) / (t30 * t35 ^ 2 + 0.1e1) * t50; 0, 0.2e1 * (-((t51 * t76 + (t56 * t47 + t74) * pkin(2)) * t48 / 0.2e1 + (t48 * t54 * t61 - t34 * t67) * pkin(1)) / t33 - (-(t40 + (-t55 * t76 - t73) * pkin(2)) * t48 / 0.2e1 + (t33 * t49 - t48 * t51) * t70) * t34 * t31) / (t31 * t34 ^ 2 + 0.1e1) * t50; 0, (0.1e1 / t45 * t76 + 0.2e1 * pkin(1) * t44 * t68) / (-t44 * t75 + 0.1e1); 0, 0;];
W = t1;
