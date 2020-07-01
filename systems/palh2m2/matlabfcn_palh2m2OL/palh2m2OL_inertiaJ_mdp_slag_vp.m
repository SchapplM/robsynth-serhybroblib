% Calculate joint inertia matrix for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh2m2OL_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 18:09
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m2OL_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(5,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'palh2m2OL_inertiaJ_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 18:06:44
% EndTime: 2020-06-30 18:06:48
% DurationCPUTime: 0.50s
% Computational Cost: add. (944->165), mult. (1774->194), div. (0->0), fcn. (2130->10), ass. (0->78)
t117 = sin(qJ(3));
t161 = cos(qJ(3));
t165 = (t161 * MDP(16) - t117 * MDP(17)) * pkin(4);
t116 = sin(qJ(4));
t160 = cos(qJ(4));
t164 = (t160 * MDP(23) - t116 * MDP(24)) * pkin(2);
t118 = sin(qJ(2));
t120 = cos(qJ(2));
t124 = t117 * t118 - t120 * t161;
t99 = t117 * t120 + t118 * t161;
t121 = t116 * t99 + t124 * t160;
t136 = -t120 * pkin(4) - pkin(1);
t90 = pkin(2) * t124 + t136;
t79 = pkin(5) * t121 + t90;
t163 = 0.2e1 * t79;
t162 = 0.2e1 * t90;
t159 = cos(qJ(5));
t158 = pkin(4) * t117;
t115 = sin(qJ(5));
t157 = t115 * pkin(5);
t87 = -t116 * t124 + t160 * t99;
t76 = t115 * t87 + t121 * t159;
t156 = t76 * MDP(34);
t155 = t76 * MDP(35);
t154 = t76 * MDP(36);
t106 = pkin(4) * t161 + pkin(2);
t95 = t160 * t106 - t116 * t158;
t92 = pkin(5) + t95;
t89 = t159 * t92;
t135 = t160 * t158;
t96 = t116 * t106 + t135;
t82 = -t115 * t96 + t89;
t153 = t82 * MDP(30);
t141 = t159 * t96;
t83 = t115 * t92 + t141;
t152 = t83 * MDP(31);
t109 = t160 * pkin(2);
t105 = t109 + pkin(5);
t100 = t159 * t105;
t146 = t116 * pkin(2);
t93 = -t115 * t146 + t100;
t151 = t93 * MDP(30);
t134 = t159 * t146;
t94 = t115 * t105 + t134;
t150 = t94 * MDP(31);
t149 = t95 * MDP(23);
t148 = t96 * MDP(24);
t114 = sin(qJ(6));
t147 = MDP(38) * t114;
t119 = cos(qJ(6));
t143 = t119 * MDP(32);
t112 = t114 ^ 2;
t113 = t119 ^ 2;
t78 = -t115 * t121 + t159 * t87;
t142 = -t76 * MDP(28) + ((t112 - t113) * MDP(33) + MDP(27)) * t78;
t140 = t114 * t119 * MDP(33);
t139 = t112 * MDP(32) + MDP(29) + 0.2e1 * t140;
t133 = MDP(22) + t139;
t132 = t87 * MDP(20) - t121 * MDP(21) + t142;
t81 = pkin(3) + t82;
t131 = t76 * t83 + t78 * t81;
t91 = pkin(3) + t93;
t130 = t76 * t94 + t78 * t91;
t129 = MDP(15) + t133;
t128 = -t114 * MDP(34) - t119 * MDP(35);
t111 = t159 * pkin(5);
t104 = t111 + pkin(3);
t127 = t104 * t78 + t157 * t76;
t126 = t99 * MDP(13) - t124 * MDP(14) + t132;
t125 = -t143 * t78 + t156;
t123 = 0.2e1 * MDP(37) * t119 - 0.2e1 * t147;
t122 = (MDP(30) * t159 - t115 * MDP(31)) * pkin(5);
t110 = pkin(3) * t119;
t102 = t104 * t119;
t88 = t91 * t119;
t80 = t81 * t119;
t72 = t76 * pkin(3) + t79;
t1 = [t121 * MDP(23) * t162 + MDP(1) + (MDP(18) * t87 - 0.2e1 * MDP(19) * t121 + MDP(24) * t162) * t87 + (MDP(11) * t99 - 0.2e1 * MDP(12) * t124) * t99 + (MDP(31) * t163 + (MDP(32) * t113 + MDP(25) - 0.2e1 * t140) * t78) * t78 + (MDP(4) * t118 + 0.2e1 * MDP(5) * t120) * t118 + 0.2e1 * (MDP(16) * t124 + t99 * MDP(17)) * t136 + 0.2e1 * (-MDP(10) * t118 + MDP(9) * t120) * pkin(1) + (MDP(30) * t163 + t72 * t123 + t154 + 0.2e1 * (-MDP(34) * t119 + MDP(35) * t114 - MDP(26)) * t78) * t76; t118 * MDP(6) + t120 * MDP(7) + (MDP(38) * t131 + t155) * t119 + (MDP(37) * t131 + t125) * t114 + t126; t81 * t123 + MDP(8) + t129 - 0.2e1 * t148 + 0.2e1 * t149 - 0.2e1 * t152 + 0.2e1 * t153 + 0.2e1 * t165; (MDP(38) * t130 + t155) * t119 + (MDP(37) * t130 + t125) * t114 + t126; (t109 + t95) * MDP(23) + (-t135 + (-pkin(2) - t106) * t116) * MDP(24) + (t100 + t89 + (-t96 - t146) * t115) * MDP(30) + (-t134 - t141 + (-t105 - t92) * t115) * MDP(31) + (t88 + t80) * MDP(37) + (-t81 - t91) * t147 + t129 + t165; t91 * t123 + t129 - 0.2e1 * t150 + 0.2e1 * t151 + 0.2e1 * t164; (MDP(38) * t127 + t155) * t119 + (MDP(37) * t127 + t125) * t114 + t132; t149 - t148 + (t111 + t82) * MDP(30) + (-t141 + (-pkin(5) - t92) * t115) * MDP(31) + (t102 + t80) * MDP(37) + (-t104 - t81) * t147 + t133; (t111 + t93) * MDP(30) + (-t134 + (-pkin(5) - t105) * t115) * MDP(31) + (t102 + t88) * MDP(37) + (-t104 - t91) * t147 + t133 + t164; t104 * t123 + 0.2e1 * t122 + t133; (pkin(3) * MDP(38) * t78 + t155) * t119 + (t156 + (pkin(3) * MDP(37) - t143) * t78) * t114 + t142; t153 - t152 + (t110 + t80) * MDP(37) + (-pkin(3) - t81) * t147 + t139; t151 - t150 + (t110 + t88) * MDP(37) + (-pkin(3) - t91) * t147 + t139; (t110 + t102) * MDP(37) + (-pkin(3) - t104) * t147 + t122 + t139; pkin(3) * t123 + t139; -t154 + (MDP(34) * t78 - MDP(37) * t72) * t119 + (-MDP(35) * t78 + MDP(38) * t72) * t114; (-MDP(38) * t83 - MDP(35)) * t119 + (-MDP(37) * t83 - MDP(34)) * t114; (-MDP(38) * t94 - MDP(35)) * t119 + (-MDP(37) * t94 - MDP(34)) * t114; (-MDP(37) * t114 - MDP(38) * t119) * t157 + t128; t128; MDP(36);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
