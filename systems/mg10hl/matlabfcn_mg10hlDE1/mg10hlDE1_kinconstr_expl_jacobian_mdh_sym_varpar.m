% Jacobian of explicit kinematic constraints of
% mg10hlDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [17x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AC,AE,CG,DC,ED,GK,GP,HP,LW,ML,OT,PM,TA,TE,phi23,phi3,phi34]';
% 
% Output:
% W [15x6]
%  Derivative of the joint coordinates w.r.t minimal coordinates
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 12:53
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function W = mg10hlDE1_kinconstr_expl_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'mg10hlDE1_kinconstr_expl_jacobian_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [17 1]), ...
  'mg10hlDE1_kinconstr_expl_jacobian_mdh_sym_varpar: pkin has to be [17x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-11 12:35:23
% EndTime: 2020-04-11 12:35:24
% DurationCPUTime: 0.47s
% Computational Cost: add. (2527->91), mult. (2317->158), div. (326->26), fcn. (870->10), ass. (0->84)
t170 = -2 * pkin(1);
t114 = sin(pkin(15));
t115 = cos(pkin(15));
t117 = sin(qJ(2));
t119 = cos(qJ(2));
t102 = t119 * t114 + t117 * t115;
t162 = t102 * pkin(1);
t151 = pkin(2) * t162;
t100 = -0.2e1 * t151;
t131 = pkin(2) ^ 2;
t132 = pkin(1) ^ 2;
t152 = t131 + t132;
t97 = t100 + t152;
t169 = 0.2e1 / t97 ^ 2;
t116 = qJ(6) + pkin(8);
t111 = 0.1e1 / t116;
t112 = 0.1e1 / t116 ^ 2;
t103 = -t117 * t114 + t119 * t115;
t168 = -0.2e1 * t103 ^ 2;
t167 = -pkin(5) - pkin(4);
t166 = -pkin(5) + pkin(4);
t137 = -qJ(6) ^ 2 + (-2 * qJ(6) - pkin(8)) * pkin(8);
t153 = (pkin(6) ^ 2 - pkin(7) ^ 2);
t105 = -t137 + t153;
t165 = -t105 / 0.2e1;
t164 = t105 / 0.4e1;
t163 = t111 / 0.2e1;
t161 = t103 * pkin(1);
t154 = t100 + t132;
t87 = (pkin(2) - t167) * (pkin(2) + t167) + t154;
t88 = (pkin(2) - t166) * (pkin(2) + t166) + t154;
t159 = t87 * t88;
t133 = sqrt(-t159);
t160 = 0.1e1 / t133 * (t87 + t88) * pkin(2) * t161;
t149 = -pkin(7) - t116;
t107 = pkin(6) - t149;
t148 = -pkin(7) + t116;
t108 = pkin(6) - t148;
t109 = pkin(6) + t148;
t110 = pkin(6) + t149;
t156 = t110 * t109;
t145 = t108 * t156;
t139 = t107 * t145;
t134 = sqrt(-t139);
t158 = 0.1e1 / t134 * (-t145 + (t156 + (t109 - t110) * t108) * t107);
t157 = t103 * t133;
t124 = 0.1e1 / pkin(7);
t155 = t124 / pkin(6);
t150 = t103 * t169;
t147 = pkin(2) * t157;
t95 = 0.1e1 / t97;
t146 = t102 * t133 * t95;
t144 = t112 * t155;
t129 = pkin(4) ^ 2;
t143 = -t129 + t152;
t138 = -t145 / 0.4e1;
t127 = pkin(5) ^ 2;
t120 = cos(pkin(16));
t118 = sin(pkin(16));
t113 = t111 * t112;
t106 = t137 + t153;
t104 = 1 / t106 ^ 2;
t99 = -pkin(2) * t102 + pkin(1);
t98 = -pkin(2) + t162;
t92 = -t127 + t129 + t97;
t91 = t100 + t127 + t143;
t90 = t127 - t143 + 0.2e1 * t151;
t89 = 0.1e1 / t90 ^ 2;
t85 = (t106 * t164 + t107 * t138) * t144;
t84 = 0.1e1 / t85 ^ 2;
t80 = (-t106 / 0.4e1 + t164) * t134 * t144;
t77 = pkin(2) * t103 * t92 + t133 * t99;
t76 = -t133 * t98 + t161 * t91;
t75 = t99 * t92 - t147;
t74 = -pkin(1) * t157 - t98 * t91;
t73 = 0.1e1 / t75 ^ 2;
t72 = 0.1e1 / t74 ^ 2;
t70 = t85 * t118 + t80 * t120;
t69 = t80 * t118 - t85 * t120;
t68 = 0.1e1 / t69 ^ 2;
t67 = (t113 * t139 / 0.2e1 + (t165 * t113 + t163) * t106 + (t116 * t165 + t138 + (t156 / 0.4e1 + (t109 / 0.4e1 - t110 / 0.4e1) * t108) * t107) * t112) * t155;
t64 = ((-t112 * t134 + t158 * t163) / t106 - (-t106 * t112 - 0.2e1 * t111 * t116) * t134 * t104) * pkin(7) * t116 * t124 / (-t104 * t139 + 0.1e1);
t62 = ((-t106 / 0.8e1 + t105 / 0.8e1) * t112 * t158 + ((t106 / 0.2e1 + t165) * t113 + t112 * t116 / 0.2e1 + t163) * t134) * t155;
t1 = [1, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0; 0, (((t131 * pkin(1) * t168 + t99 * t160) * t95 + ((-t102 * t92 - t157) * t95 + t77 * pkin(1) * t150) * pkin(2)) / t75 - (t146 + ((t170 * t99 - t160 - t92) * t95 + t75 * pkin(1) * t169) * t103) * pkin(2) * t77 * t73) / (t77 ^ 2 * t73 + 0.1e1) * t97, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t64; 0, 0, 0, 0, 0, ((t67 * t118 + t62 * t120) / t69 - (t62 * t118 - t67 * t120) * t70 * t68) / (t70 ^ 2 * t68 + 0.1e1); 0, 0, 1, 0, 0, 0; 0, 0, 0, 1, 0, 0; 0, 0, 0, 0, 1, 0; 0, 0, 0, 0, 0, t64; 0, (((t132 * pkin(2) * t168 - t98 * t160) * t95 + ((-t102 * t91 - t157) * t95 + t76 * pkin(2) * t150) * pkin(1)) / t74 - (t146 + ((0.2e1 * t98 * pkin(2) - t160 - t91) * t95 + t74 * pkin(2) * t169) * t103) * pkin(1) * t76 * t72) / (t76 ^ 2 * t72 + 0.1e1) * t97, 0, 0, 0, 0; 0, (0.1e1 / t90 * t160 + t89 * t147 * t170) / (-t159 * t89 + 0.1e1), 0, 0, 0, 0; 0, 0, 0, 0, 0, (t62 / t85 - t67 * t80 * t84) / (t80 ^ 2 * t84 + 0.1e1); 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
W = t1;
