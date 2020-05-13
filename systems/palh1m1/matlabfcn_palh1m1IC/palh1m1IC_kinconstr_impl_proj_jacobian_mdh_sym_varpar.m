% Jacobian of implicit kinematic constraints of
% palh1m1IC
% projection from active to passive joints coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% 
% Output:
% B21 [(no of passive joints)x(no of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 20:03
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function B21 = palh1m1IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: pkin has to be [20x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_projection_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:56:38
% EndTime: 2020-04-15 19:56:38
% DurationCPUTime: 0.15s
% Computational Cost: add. (445->39), mult. (370->66), div. (30->7), fcn. (260->20), ass. (0->43)
t124 = -qJ(7) + pkin(19);
t119 = -qJ(10) + t124;
t111 = sin(t119);
t112 = cos(t119);
t120 = qJ(3) + qJ(4) + pkin(18);
t113 = sin(t120);
t115 = cos(t120);
t141 = (t111 * t113 - t112 * t115) * pkin(10);
t121 = pkin(20) + qJ(7) + qJ(2);
t114 = sin(t121);
t137 = pkin(3) * t114;
t103 = t137 + sin(qJ(2)) * pkin(1);
t116 = cos(t121);
t136 = pkin(3) * t116;
t104 = t136 + cos(qJ(2)) * pkin(1);
t128 = sin(qJ(6));
t130 = cos(qJ(6));
t97 = 0.1e1 / (-t114 * t128 - t116 * t130) / pkin(7) / pkin(3);
t140 = (t103 * t128 + t104 * t130) * t97 * pkin(7);
t126 = qJ(8) + qJ(9);
t122 = sin(t126);
t127 = sin(qJ(8));
t105 = pkin(2) * t127 - pkin(12) * t122;
t123 = cos(t126);
t129 = cos(qJ(8));
t106 = -pkin(2) * t129 + pkin(12) * t123;
t139 = pkin(6) / (t105 * t123 + t106 * t122) / pkin(12);
t109 = pkin(10) * t113;
t110 = pkin(10) * t115;
t108 = t112 * pkin(8);
t100 = t108 - pkin(4) * cos(t124);
t107 = t111 * pkin(8);
t99 = t107 - pkin(4) * sin(t124);
t135 = t100 * t110 - t99 * t109;
t134 = -t100 * t111 + t112 * t99;
t96 = 0.1e1 / pkin(8) / t141;
t131 = t96 * t140;
t125 = qJ(3) + pkin(17);
t118 = cos(t125);
t117 = sin(t125);
t102 = t110 + cos(qJ(3)) * pkin(5);
t101 = t109 + sin(qJ(3)) * pkin(5);
t1 = [0, t134 * pkin(8) * t131, -(t101 * t111 - t102 * t112) * t96 * pkin(8), 0; 0, (t103 * t116 - t104 * t114) * t97 * pkin(3), 0, 0; 0, t140, 0, 0; 0, 0, (-t117 * t122 - t118 * t123) * pkin(12) * t139, 0; 0, 0, (-t105 * t117 + t106 * t118) * t139, 0; 0, t135 * t131, -(-t101 * t115 + t102 * t113) * t96 * pkin(10), 0; 0, -0.1e1 + ((-pkin(7) * t130 - t137) * t104 - (pkin(7) * t128 - t136) * t103) * t97, 0, 0; 0, 0, -0.1e1 + (-t117 * t127 - t118 * t129) * pkin(2) * t139, 0; 0, ((-t134 + t141) * pkin(8) + t135) * t131, -0.1e1 - ((t109 + t108) * t102 - (t110 + t107) * t101) * t96, 0;];
B21 = t1;
