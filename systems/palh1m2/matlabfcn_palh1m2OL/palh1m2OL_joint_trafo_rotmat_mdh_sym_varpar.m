% Calculate homogenous joint transformation matrices for
% palh1m2OL
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
% T_mdh [4x4x16]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(16+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-29 20:34
% Revision: bc59515823ab4a8d0fec19bf3bf92c32c39a66b0 (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = palh1m2OL_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [20x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-29 20:25:06
% EndTime: 2020-06-29 20:25:06
% DurationCPUTime: 0.13s
% Computational Cost: add. (33->29), mult. (40->24), div. (0->0), fcn. (108->34), ass. (0->43)
t130 = cos(qJ(1));
t129 = cos(qJ(2));
t128 = cos(qJ(3));
t127 = cos(qJ(4));
t126 = cos(qJ(5));
t125 = cos(qJ(6));
t124 = cos(qJ(7));
t123 = cos(qJ(8));
t122 = cos(qJ(9));
t121 = sin(qJ(1));
t120 = sin(qJ(2));
t119 = sin(qJ(3));
t118 = sin(qJ(4));
t117 = sin(qJ(5));
t116 = sin(qJ(6));
t115 = sin(qJ(7));
t114 = sin(qJ(8));
t113 = sin(qJ(9));
t112 = cos(qJ(10));
t111 = cos(qJ(11));
t110 = cos(qJ(12));
t109 = cos(qJ(13));
t108 = sin(qJ(10));
t107 = sin(qJ(11));
t106 = sin(qJ(12));
t105 = sin(qJ(13));
t104 = cos(pkin(17));
t103 = cos(pkin(18));
t102 = cos(pkin(19));
t101 = cos(pkin(20));
t100 = sin(pkin(17));
t99 = sin(pkin(18));
t98 = sin(pkin(19));
t97 = sin(pkin(20));
t96 = t100 * t106 - t104 * t110;
t95 = t103 * t109 - t99 * t105;
t94 = -t102 * t112 - t98 * t108;
t93 = -t101 * t111 + t97 * t107;
t92 = -t100 * t110 - t104 * t106;
t91 = t103 * t105 + t99 * t109;
t90 = -t102 * t108 + t98 * t112;
t89 = -t101 * t107 - t97 * t111;
t1 = [t130, -t121, 0, 0; t121, t130, 0, 0; 0, 0, 1, pkin(13); -t120, -t129, 0, pkin(15); 0, 0, -1, 0; t129, -t120, 0, 0; t119, t128, 0, pkin(1); -t128, t119, 0, 0; 0, 0, 1, 0; t127, -t118, 0, pkin(5); t118, t127, 0, 0; 0, 0, 1, 0; t126, -t117, 0, pkin(9); 0, 0, -1, -pkin(11); t117, t126, 0, 0; t125, -t116, 0, -pkin(14); 0, 0, -1, 0; t116, t125, 0, -pkin(16); t124, -t115, 0, pkin(1); t115, t124, 0, 0; 0, 0, 1, 0; t123, -t114, 0, 0; t114, t123, 0, 0; 0, 0, 1, 0; -t122, t113, 0, pkin(2); -t113, -t122, 0, 0; 0, 0, 1, 0; t94, -t90, 0, t102 * pkin(4); t90, t94, 0, -t98 * pkin(4); 0, 0, 1, 0; t93, -t89, 0, t101 * pkin(3); t89, t93, 0, t97 * pkin(3); 0, 0, 1, 0; t96, -t92, 0, t104 * pkin(6); t92, t96, 0, t100 * pkin(6); 0, 0, 1, 0; t95, -t91, 0, t103 * pkin(10); t91, t95, 0, t99 * pkin(10); 0, 0, 1, 0; 1, 0, 0, pkin(7); 0, 1, 0, 0; 0, 0, 1, 0; 1, 0, 0, pkin(12); 0, 1, 0, 0; 0, 0, 1, 0; -1, 0, 0, pkin(8); 0, -1, 0, 0; 0, 0, 1, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,16);             % numerisch
else,                         T_mdh = sym('xx', [4,4,16]); end % symbolisch

for i = 1:16
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
