% Jacobian of explicit kinematic constraints of
% palh4m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AD,CB,CE,EP,HC,OT,TA,TD]';
% 
% Output:
% W [9x5]
%  Derivative of the joint coordinates w.r.t minimal coordinates
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 21:48
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function W = palh4m1TE_kinconstr_expl_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh4m1TE_kinconstr_expl_jacobian_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'palh4m1TE_kinconstr_expl_jacobian_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-11 21:35:19
% EndTime: 2020-04-11 21:35:20
% DurationCPUTime: 0.27s
% Computational Cost: add. (1533->59), mult. (1558->126), div. (68->14), fcn. (406->4), ass. (0->62)
t98 = -2 * pkin(2);
t97 = 2 * pkin(2);
t61 = (qJ(2) + pkin(6));
t96 = -2 * t61;
t63 = cos(qJ(5));
t95 = pkin(2) * t63;
t62 = sin(qJ(5));
t93 = t62 * pkin(1);
t92 = t63 * pkin(1);
t82 = pkin(2) * t93;
t56 = -0.2e1 * t82;
t66 = pkin(2) ^ 2;
t67 = pkin(1) ^ 2;
t83 = t66 + t67;
t53 = t56 + t83;
t57 = t61 ^ 2;
t64 = pkin(3) ^ 2;
t50 = -t57 + t64 + t53;
t55 = -pkin(2) * t62 + pkin(1);
t79 = -pkin(3) - t61;
t84 = t56 + t67;
t44 = ((pkin(2) - t79) * (pkin(2) + t79)) + t84;
t78 = pkin(3) - t61;
t45 = ((pkin(2) - t78) * (pkin(2) + t78)) + t84;
t88 = t44 * t45;
t68 = sqrt(-t88);
t85 = t63 * t68;
t81 = pkin(2) * t85;
t36 = t55 * t50 - t81;
t34 = 0.1e1 / t36 ^ 2;
t38 = t50 * t95 + t55 * t68;
t91 = t34 * t38;
t43 = 0.1e1 / t68;
t90 = 0.2e1 * t43 * (-t44 * t78 - t79 * t45);
t89 = t43 * (t45 + t44) * pkin(2) * t92;
t87 = t61 * t63;
t86 = t62 * t68;
t52 = 0.1e1 / t53 ^ 2;
t80 = t52 * t95;
t76 = -t64 + t83;
t49 = t56 + t57 + t76;
t54 = -pkin(2) + t93;
t35 = -pkin(1) * t85 - t54 * t49;
t32 = 0.1e1 / t35 ^ 2;
t37 = t49 * t92 - t54 * t68;
t77 = 0.1e1 / (t37 ^ 2 * t32 + 0.1e1) * t53 * t61;
t75 = -t90 / 0.2e1;
t74 = t90 / 0.2e1;
t71 = 0.1e1 / t35 * t77;
t70 = t63 * t75;
t69 = t32 * t37 * t77;
t60 = t63 ^ 2;
t59 = 1 / t61 ^ 2;
t58 = 1 / t61;
t51 = 0.1e1 / t53;
t48 = t57 - t76 + 0.2e1 * t82;
t47 = 0.1e1 / t48 ^ 2;
t46 = 0.1e1 / t48;
t39 = 0.1e1 / (-t47 * t88 + 0.1e1);
t33 = 0.1e1 / t36;
t30 = 0.1e1 / (t38 ^ 2 * t34 + 0.1e1);
t1 = [1, 0, 0, 0, 0; 0, (((0.2e1 * pkin(1) * t87 + t54 * t75) * t58 - t37 * t59) * t71 - ((pkin(1) * t70 + t54 * t96) * t58 - t35 * t59) * t69) * t51, 0, 0, (((t67 * t60 * t98 - t54 * t89) * t51 + ((-t62 * t49 - t85) * t51 + 0.2e1 * t37 * t80) * pkin(1)) * t71 - (t51 * t86 + ((t54 * t97 - t49 - t89) * t51 + t35 * t52 * t97) * t63) * pkin(1) * t69) * t58; 0, 1, 0, 0, 0; 0, ((t58 * t74 - t59 * t68) * t46 - (-t48 * t59 + 0.2e1) * t68 * t47) * t39 * t61, 0, 0, (-0.2e1 * pkin(1) * t47 * t81 + t46 * t89) * t39; 0, 0, 1, 0, 0; 0, 0, 0, 1, 0; 0, 0, 0, 0, 1; 0, (-(t55 * t74 + t87 * t98) * t33 + (pkin(2) * t70 + t55 * t96) * t91) * t30, 0, 0, 0.2e1 * (-((t55 * t89 + (-t62 * t50 - t85) * pkin(2)) * t51 / 0.2e1 + (-t66 * t60 * t51 + t38 * t80) * pkin(1)) * t33 - (-(t86 + (-t50 - t89) * t63) * t51 / 0.2e1 + (-t36 * t52 + t51 * t55) * t92) * pkin(2) * t91) * t30 * t53; 0, 0, 0, 0, 0;];
W = t1;
