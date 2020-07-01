% Calculate Gravitation load on the joints for
% fourbar1turnDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar1turnDE1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:36
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1turnDE1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'fourbar1turnDE1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:35:29
% EndTime: 2020-06-27 16:35:36
% DurationCPUTime: 1.30s
% Computational Cost: add. (8844->96), mult. (12864->184), div. (612->9), fcn. (3605->10), ass. (0->86)
t155 = cos(qJ(2));
t153 = sin(qJ(2));
t163 = pkin(2) ^ 2;
t164 = pkin(1) ^ 2;
t208 = pkin(2) * t155;
t194 = -0.2e1 * pkin(1) * t208 + t164;
t148 = t163 + t194;
t222 = pkin(3) ^ 2;
t223 = pkin(4) ^ 2;
t193 = t222 - t223;
t143 = t148 + t193;
t140 = pkin(1) * t153 * t143;
t150 = pkin(1) * t155 - pkin(2);
t217 = -pkin(3) - pkin(4);
t141 = (pkin(2) - t217) * (pkin(2) + t217) + t194;
t216 = pkin(4) - pkin(3);
t142 = (pkin(2) - t216) * (pkin(2) + t216) + t194;
t167 = sqrt(-t141 * t142);
t200 = t155 * t167;
t209 = pkin(2) * t153;
t191 = pkin(1) * t209;
t203 = 0.1e1 / t167 * (-t141 - t142) * t191;
t118 = t140 + (-t200 + (-0.2e1 * t150 * pkin(2) - t203) * t153) * pkin(1);
t135 = -t150 * t167 + t140;
t198 = t118 - t135;
t185 = t198 * t153;
t201 = t153 * t167;
t152 = t153 ^ 2;
t218 = 0.2e1 * t152;
t120 = -t150 * t203 + t164 * pkin(2) * t218 + (t143 * t155 + t201) * pkin(1);
t132 = -pkin(1) * t201 - t143 * t150;
t197 = t120 + t132;
t224 = t197 * t155 + t185;
t154 = sin(qJ(1));
t213 = g(2) * t154;
t156 = cos(qJ(1));
t214 = g(1) * t156;
t183 = t213 + t214;
t176 = 0.2e1 * t183;
t146 = 0.1e1 / t148 ^ 2;
t159 = 0.1e1 / t223;
t144 = t148 - t193;
t149 = pkin(1) - t208;
t133 = -pkin(2) * t201 + t144 * t149;
t139 = t144 * t209;
t134 = t149 * t167 + t139;
t195 = t133 ^ 2 + t134 ^ 2;
t126 = t195 * t159 * t146;
t122 = t126 ^ (-0.1e1 / 0.2e1);
t145 = 0.1e1 / t148;
t221 = t122 * t145;
t162 = 0.1e1 / t222;
t196 = t132 ^ 2 + t135 ^ 2;
t127 = t196 * t162 * t146;
t124 = t127 ^ (-0.1e1 / 0.2e1);
t220 = t124 * t145;
t186 = t146 * t191;
t215 = g(1) * t154;
t212 = g(2) * t156;
t211 = g(3) * t133;
t210 = g(3) * t134;
t207 = t132 * t153;
t206 = t132 * t155;
t205 = t135 * t153;
t204 = t135 * t155;
t202 = t153 * t155;
t161 = 0.1e1 / pkin(3);
t187 = t161 * t220;
t184 = t155 * t187;
t199 = (t135 * t184 + t187 * t207) * t156;
t192 = pkin(1) * pkin(2) * t146;
t190 = 0.2e1 * t146;
t119 = t139 + (-t200 + (0.2e1 * t149 * pkin(1) - t203) * t153) * pkin(2);
t121 = t149 * t203 + t163 * pkin(1) * t218 + (t144 * t155 + t201) * pkin(2);
t181 = 0.4e1 * t145 * t186;
t189 = ((t119 * t133 + t121 * t134) * t190 - t195 * t181) * t159 / t126 * t221;
t188 = ((t118 * t132 + t120 * t135) * t190 - t196 * t181) * t162 / t127 * t220;
t180 = -t214 / 0.2e1 - t213 / 0.2e1;
t179 = t132 * t152 + t135 * t202;
t178 = t206 / 0.2e1 - t205 / 0.2e1;
t177 = -t204 / 0.2e1 - t207 / 0.2e1;
t173 = 0.2e1 * t132 * t202 - 0.2e1 * t135 * t152;
t172 = t197 * t153 - t198 * t155;
t158 = 0.1e1 / pkin(4);
t115 = t154 * t132 * t184;
t1 = [t183 * MDP(3) + (-g(1) * t115 + (t205 * t215 - (t205 - t206) * t212) * t187) * MDP(17) + (-g(2) * t199 - (-t204 - t207) * t187 * t215) * MDP(18) + (-(t133 * MDP(24) + t134 * MDP(25)) * t158 * t221 - t153 * MDP(10) + t155 * MDP(9) + MDP(2)) * (-t212 + t215); (-g(3) * t155 + t183 * t153) * MDP(9) + (g(3) * t153 + t183 * t155) * MDP(10) + (-g(1) * t199 + ((g(3) * t177 - t183 * t178) * t188 + ((-0.2e1 * g(3) * t179 - t183 * t173) * t192 + (g(3) * t224 - (-t118 * t155 + t120 * t153) * t214 - t172 * t213) * t145) * t124) * t161) * MDP(17) + (-g(2) * t115 + ((-g(3) * t178 - t183 * t177) * t188 + ((-g(3) * t173 + t179 * t176) * t192 + (-g(3) * t172 - (t120 * t155 + t185) * t213 - t224 * t214) * t145) * t124) * t161) * MDP(18) + (((t210 / 0.2e1 + t180 * t133) * t189 + ((-g(3) * t121 + t183 * t119) * t145 + (-t133 * t176 + 0.2e1 * t210) * t186) * t122) * MDP(24) + ((-t211 / 0.2e1 + t180 * t134) * t189 + ((g(3) * t119 + t183 * t121) * t145 + (-t134 * t176 - 0.2e1 * t211) * t186) * t122) * MDP(25)) * t158;];
taug = t1;
