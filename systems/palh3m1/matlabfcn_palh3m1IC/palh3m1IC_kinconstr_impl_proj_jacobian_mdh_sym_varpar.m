% Jacobian of implicit kinematic constraints of
% palh3m1IC
% projection from active to passive joints coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% 
% Output:
% B21 [(no of passive joints)x(no of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 17:32
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function B21 = palh3m1IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: qJ has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_projection_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:27:33
% EndTime: 2020-04-20 17:27:34
% DurationCPUTime: 0.13s
% Computational Cost: add. (387->31), mult. (302->48), div. (24->4), fcn. (210->14), ass. (0->32)
t91 = -qJ(7) + pkin(15);
t89 = -qJ(8) + t91;
t86 = cos(t89);
t100 = t86 * pkin(7);
t84 = sin(t89);
t101 = t84 * pkin(7);
t88 = qJ(3) + qJ(4) + pkin(14);
t83 = cos(t88);
t102 = pkin(9) * t83;
t82 = sin(t88);
t103 = pkin(9) * t82;
t99 = t100 * t103 + t101 * t102;
t71 = 0.1e1 / t99;
t90 = pkin(16) + qJ(7) + qJ(2);
t85 = sin(t90);
t87 = cos(t90);
t92 = sin(qJ(6));
t93 = cos(qJ(6));
t72 = 0.1e1 / (-t85 * t93 + t87 * t92) / pkin(5) / pkin(2);
t106 = pkin(2) * t85;
t80 = -t106 - sin(qJ(2)) * pkin(1);
t105 = pkin(2) * t87;
t81 = t105 + cos(qJ(2)) * pkin(1);
t95 = (t93 * t80 + t92 * t81) * t72 * pkin(5);
t94 = t71 * t95;
t74 = -t101 + pkin(3) * sin(t91);
t75 = -t100 + pkin(3) * cos(t91);
t107 = (t74 * t83 + t75 * t82) * pkin(9);
t98 = -t100 * t74 + t75 * t101;
t79 = t102 + cos(qJ(3)) * pkin(4);
t78 = -t103 - sin(qJ(3)) * pkin(4);
t1 = [0, t98 * t94, (t86 * t78 - t84 * t79) * t71 * pkin(7), 0; 0, (-t87 * t80 - t85 * t81) * t72 * pkin(2), 0, 0; 0, -t95, 0, 0; 0, -t94 * t107, (t83 * t78 + t82 * t79) * t71 * pkin(9), 0; 0, -0.1e1 + ((t93 * pkin(5) - t105) * t80 - (-t92 * pkin(5) + t106) * t81) * t72, 0, 0; 0, -(t98 + t99 + t107) * t94, -0.1e1 + (-(t100 - t102) * t78 + (t101 + t103) * t79) * t71, 0;];
B21 = t1;
