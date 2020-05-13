% Jacobian of implicit kinematic constraints of
% palh1m2IC
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
% Datum: 2020-05-02 23:49
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function B21 = palh1m2IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: pkin has to be [20x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_projection_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:42:15
% EndTime: 2020-05-02 23:42:15
% DurationCPUTime: 0.09s
% Computational Cost: add. (230->46), mult. (57->46), div. (33->10), fcn. (53->25), ass. (0->38)
t139 = pkin(6) / sin(qJ(9));
t133 = qJ(4) + pkin(18);
t134 = (pkin(19) + qJ(3));
t125 = (t133 + t134);
t136 = (qJ(6) - qJ(2));
t123 = (t125 + t136);
t118 = -2 * qJ(7) - pkin(20) + t123;
t124 = (t125 - t136);
t119 = pkin(20) + t124;
t138 = cos((qJ(10) - t118)) + cos((qJ(10) - t119));
t116 = 0.1e1 / pkin(3);
t137 = t116 / t138;
t135 = qJ(10) + qJ(7);
t109 = pkin(17) + qJ(3);
t132 = qJ(7) + pkin(20);
t126 = t132 - t136;
t105 = cos(t126);
t102 = 0.1e1 / t105;
t131 = pkin(1) * t102 / pkin(7);
t101 = cos((t125 - t135));
t100 = 0.1e1 / t101;
t114 = 0.1e1 / pkin(8);
t130 = pkin(5) * t100 * t114;
t129 = 0.1e1 / pkin(2) * t139;
t128 = -qJ(8) + t109;
t113 = 0.1e1 / pkin(10);
t127 = pkin(4) * t113 * t137;
t120 = -qJ(7) + t123;
t121 = -qJ(7) + t124;
t122 = (-cos(t118) - cos(t119)) * pkin(3) + (-cos(t120) - cos(t121)) * pkin(1);
t112 = 0.1e1 / pkin(12);
t108 = cos(t136);
t107 = sin(t133);
t106 = sin(t132);
t104 = cos(qJ(9) - t128);
t103 = cos((t134 - t135));
t92 = pkin(1) * (sin((qJ(10) - t136)) + sin((qJ(10) + t136))) + (sin((qJ(10) - t126)) + sin((qJ(10) + t126))) * pkin(3);
t1 = [0, -t92 * t127, (-pkin(5) * t103 - pkin(10) * t101) * t113 * t100, 0; 0, t106 * t131, 0, 0; 0, (-pkin(1) * t108 - pkin(3) * t105) * t116 * t102, 0, 0; 0, 0, t104 * t129, 0; 0, 0, (-pkin(12) * t104 + pkin(2) * cos(t128)) * t112 * t129, 0; 0, (t122 * pkin(4) + (t138 * pkin(3) + (cos((-qJ(10) + t120)) + cos((-qJ(10) + t121))) * pkin(1)) * pkin(8)) * t114 * t137, t107 * t130, 0; 0, (pkin(3) * t106 + pkin(7) * t108) * t116 * t131, 0, 0; 0, 0, -0.1e1 + (cos(t109) * cos(qJ(8)) + sin(t109) * sin(qJ(8))) * t112 * t139, 0; 0, (t92 * pkin(8) + t122 * pkin(10)) * t114 * t127, (pkin(8) * t103 + pkin(10) * t107) * t113 * t130, 0;];
B21 = t1;
